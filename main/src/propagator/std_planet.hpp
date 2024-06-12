
/*! @file
 * @brief A Propagator class for Protoplanetary disk
 *
 * @author Noah Kubli
 */

#pragma once

#include <variant>

#include "cstone/fields/field_get.hpp"
#include "sph/particles_data.hpp"
#include "sph/sph.hpp"
#include "exchangeStarPosition.hpp"
#include "computeCentralForce.hpp"
#include "star_data.hpp"
#include "accretion.hpp"
#include "betaCooling.hpp"

#include "ipropagator.hpp"
#include "gravity_wrapper.hpp"

namespace sphexa
{

using namespace sph;
using util::FieldList;

template<class DomainType, class DataType>
class PlanetProp : public HydroProp<DomainType, DataType>
{
protected:
    using Base = HydroProp<DomainType, DataType>; // Propagator<DomainType, DataType>;
    using Base::timer;

    using T             = typename DataType::RealType;
    using KeyType       = typename DataType::KeyType;
    using Tmass         = typename DataType::HydroData::Tmass;
    using MultipoleType = ryoanji::CartesianQuadrupole<Tmass>;

    using Acc       = typename DataType::AcceleratorType;
    using MHolder_t = typename cstone::AccelSwitchType<Acc, MultipoleHolderCpu, MultipoleHolderGpu>::template type<
        MultipoleType, DomainType, typename DataType::HydroData>;

    MHolder_t mHolder_;

    StarData star;
    // std::array<double, 3> stellar_position;
    // std::array<double, 3> stellar_position_m1;

    /*! @brief the list of conserved particles fields with values preserved between iterations
     *
     * x, y, z, h and m are automatically considered conserved and must not be specified in this list
     */
    using ConservedFields = FieldList<"u", "vx", "vy", "vz", "x_m1", "y_m1", "z_m1", "du_m1", "ParticleIDs">;

    //! @brief the list of dependent particle fields, these may be used as scratch space during domain sync
    using DependentFields =
        FieldList<"rho", "p", "c", "ax", "ay", "az", "du", "keyTypeSwap", "c11", "c12", "c13", "c22", "c23", "c33", "nc">;

public:
    PlanetProp(std::ostream& output, size_t rank, const InitSettings& settings)
        : Base(output, rank)
    {
    }

    std::vector<std::string> conservedFields() const override
    {
        std::vector<std::string> ret{"x", "y", "z", "h", "m"};
        for_each_tuple([&ret](auto f) { ret.push_back(f.value); }, make_tuple(ConservedFields{}));
        return ret;
    }
    void load(const std::string& initCond, IFileReader* reader) override
    {
        // Read star position from hdf5 File
        std::string path = removeModifiers(initCond);
        if (std::filesystem::exists(path))
        {
            int snapshotIndex = numberAfterSign(initCond, ":");
            reader->setStep(path, snapshotIndex, FileMode::independent);
            star.loadOrStoreAttributes(reader);
            reader->closeStep();
            printf("star position: %lf\t%lf\t%lf\n", star.position[0], star.position[1], star.position[2]);
            printf("star mass: %lf\n", star.m);
        }
    }
    void save(IFileWriter* writer) override { star.loadOrStoreAttributes(writer); }

    void activateFields(DataType& simData) override
    {
        auto& d = simData.hydro;

        //! @brief Fields accessed in domain sync are not part of extensible lists.
        d.setConserved("x", "y", "z", "h", "m");
        d.setDependent("keys");
        std::apply([&d](auto... f) { d.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([&d](auto... f) { d.setDependent(f.value...); }, make_tuple(DependentFields{}));

        d.devData.setConserved("x", "y", "z", "h", "m");
        d.devData.setDependent("keys");
        std::apply([&d](auto... f) { d.devData.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([&d](auto... f) { d.devData.setDependent(f.value...); }, make_tuple(DependentFields{}));
    }

    void sync(DomainType& domain, DataType& simData) override
    {
        auto& d = simData.hydro;
        if (d.g != 0.0)
        {
            domain.syncGrav(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d), get<"m">(d),
                            get<ConservedFields>(d), get<DependentFields>(d));
        }
        else
        {
            domain.sync(get<"keys">(d), get<"x">(d), get<"y">(d), get<"z">(d), get<"h">(d),
                        std::tuple_cat(std::tie(get<"m">(d)), get<ConservedFields>(d)), get<DependentFields>(d));
        }
        d.treeView = domain.octreeProperties();
    }

    void computeForces(DomainType& domain, DataType& simData)
    {
        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();
        auto&  d     = simData.hydro;

        resizeNeighbors(d, domain.nParticles() * d.ngmax);
        findNeighborsSfc(first, last, d, domain.box());
        timer.step("FindNeighbors");

        computeDensity(first, last, d, domain.box());
        timer.step("Density");
        computeEOS_HydroStd(first, last, d);
        timer.step("EquationOfState");

        domain.exchangeHalos(get<"vx", "vy", "vz", "rho", "p", "c">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeIAD(first, last, d, domain.box());
        timer.step("IAD");

        domain.exchangeHalos(get<"c11", "c12", "c13", "c22", "c23", "c33">(d), get<"ax">(d), get<"ay">(d));
        timer.step("mpi::synchronizeHalos");

        computeMomentumEnergySTD(first, last, d, domain.box());
        timer.step("MomentumEnergyIAD");

        if (d.g != 0.0)
        {
            mHolder_.upsweep(d, domain);
            timer.step("Upsweep");
            //mHolder_.traverse(d, domain);
            //timer.step("Gravity");
        }
    }

    void step(DomainType& domain, DataType& simData) override
    {
        timer.start();
        auto& d = simData.hydro;

        using KeyType = typename DataType::HydroData::KeyType;

        size_t first = domain.startIndex();
        size_t last  = domain.endIndex();

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        fill(get<"keys">(d), first, last, KeyType{0});

        //planet::computeAccretionCondition(first, last, d, star);

       // planet::computeNewOrder(first, last, d, star);
        //planet::applyNewOrder<ConservedFields, DependentFields>(first, last, d, star);

        //planet::sumAccretedMassAndMomentum<DependentFields>(first, last, d, star);
        //planet::exchangeAndAccreteOnStar(star, d.minDt_m1, rank);

        domain.setEndIndex(last - star.n_accreted_local - star.n_removed_local);

        timer.step("accreteParticles");

        sync(domain, simData);
        timer.step("domain::sync");

        d.resize(domain.nParticlesWithHalos());
        domain.exchangeHalos(std::tie(get<"m">(d)), get<"ax">(d), get<"ay">(d));
        first = domain.startIndex();
        last  = domain.endIndex();

        computeForces(domain, simData);

        //planet::betaCooling(d, first, last, star);
        timer.step("betaCooling");

        planet::computeCentralForce(simData.hydro, first, last, star);
        timer.step("computeCentralForce");
        double dtu = planet::computeHeatingTimestep(d, first, last);
        computeTimestep(first, last, d, dtu);

        timer.step("Timestep");

        computePositions(first, last, d, domain.box());
        updateSmoothingLength(first, last, d);
        timer.step("UpdateQuantities");

        planet::computeAndExchangeStarPosition(star, d.minDt, d.minDt_m1, rank);
        timer.step("computeAndExchangeStarPosition");

        if (rank == 0)
        {
            printf("star position: %lf\t%lf\t%lf\n", star.position[0], star.position[1], star.position[2]);
            printf("star mass: %lf\n", star.m);
            printf("additional pot. erg.: %lf\n", star.potential);
            printf("heating timestep: %lf\t courant: %lf\n", dtu, d.minDtCourant);
        }
        printf("rank: %d, accreted %zu, removed %zu\n", rank, star.n_accreted_local, star.n_removed_local);
    }

    void saveFields(IFileWriter* writer, size_t first, size_t last, DataType& simData,
                    const cstone::Box<T>& /*box*/) override
    {
        auto output = [&](auto& d)
        {
            auto fieldPointers = d.data();
            auto indicesDone   = d.outputFieldIndices;
            auto namesDone     = d.outputFieldNames;

            for (int i = int(indicesDone.size()) - 1; i >= 0; --i)
            {
                int fidx = indicesDone[i];
                if (d.isAllocated(fidx))
                {
                    int column = std::find(d.outputFieldIndices.begin(), d.outputFieldIndices.end(), fidx) -
                                 d.outputFieldIndices.begin();
                    transferToHost(d, first, last, {d.fieldNames[fidx]});
                    std::visit([writer, c = column, key = namesDone[i]](auto field)
                               { writer->writeField(key, field->data(), c); },
                               fieldPointers[fidx]);
                    indicesDone.erase(indicesDone.begin() + i);
                    namesDone.erase(namesDone.begin() + i);
                }
            }

            if (!indicesDone.empty() && Base::rank_ == 0)
            {
                std::cout << "WARNING: the following fields are not in use and therefore not output: ";
                for (int fidx = 0; fidx < indicesDone.size() - 1; ++fidx)
                {
                    std::cout << d.fieldNames[fidx] << ",";
                }
                std::cout << d.fieldNames[indicesDone.back()] << std::endl;
            }
        };

        output(simData.hydro);
        output(simData.chem);
    }
};

} // namespace sphexa