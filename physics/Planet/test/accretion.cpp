//
// Created by Noah Kubli on 24.03.2024.
//

#include "accretion.hpp"
#include "cstone/domain/domain.hpp"
#include "cstone/fields/field_get.hpp"
#include "gtest/gtest.h"
#include "sph/particles_data.hpp"
#include "star_data.hpp"

#include <memory>
#include <mpi.h>

static void fill_grid(auto& x, auto& y, auto& z, double start, double end, size_t n_elements)
{
    assert(x.size() == y.size() == z.size());

    for (size_t i = 0; i < x.size(); i++)
    {
        const size_t iz = i / (n_elements * n_elements);
        const size_t iy = i / n_elements % n_elements;
        const size_t ix = i % n_elements;
        x[i]            = start + ix * (end - start) / n_elements;
        y[i]            = start + iy * (end - start) / n_elements;
        z[i]            = start + iz * (end - start) / n_elements;
    }
}

struct AccretionTest : public ::testing::Test
{
    using KeyType         = uint64_t;
    using T               = double;
    using ConservedFields = util::FieldList<"vx", "vy", "vz", "m">;
    using DependentFields = util::FieldList<"keys", "u", "rho">;

    inline constexpr static size_t n_elements  = 20;
    inline constexpr static size_t n_elements3 = n_elements * n_elements * n_elements;
    inline constexpr static double inner_limit = 0.2;
    inline constexpr static double dt          = 1.;

    sphexa::ParticlesData<cstone::CpuTag>       data;
    std::unique_ptr<cstone::Domain<KeyType, T>> domain_ptr = nullptr;

    StarData star;
    double   momentum_inside_x = 0.;
    double   momentum_inside_y = 0.;
    double   momentum_inside_z = 0.;
    double   mass_inside       = 0.;
    size_t   n_inside          = 0;
    int      rank = 0, numRanks = 0;

    AccretionTest()
    {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (initialized == 0) { MPI_Init(0, NULL); }
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

        int   bucketSize = 1;
        float theta      = 1.0;

        domain_ptr = std::make_unique<cstone::Domain<KeyType, T>>(rank, numRanks, bucketSize, bucketSize, theta);

        star.inner_size = inner_limit;
    }



    void initData()
    {
        data.setConserved("x", "y", "z", "h");
        std::apply([this](auto... f) { data.setConserved(f.value...); }, make_tuple(ConservedFields{}));
        std::apply([this](auto... f) { data.setDependent(f.value...); }, make_tuple(DependentFields{}));

        data.resize(n_elements3);
        std::cout << "size: " << data.x.size() << std::endl;

        fill_grid(data.x, data.y, data.z, -1., 1., n_elements);
        fill_grid(data.vx, data.vy, data.vz, 0., 1., n_elements);

        std::fill(data.h.begin(), data.h.end(), 0.05);
        std::fill(data.m.begin(), data.m.end(), 1.);

        momentum_inside_x = 0.;
        momentum_inside_y = 0.;
        momentum_inside_z = 0.;
        mass_inside       = 0.;
        n_inside          = 0;

        for (size_t i = 0; i < data.x.size(); i++)
        {
            T x    = data.x[i];
            T y    = data.y[i];
            T z    = data.z[i];
            T dist = std::sqrt(x * x + y * y + z * z);

            if (dist < inner_limit)
            {
                momentum_inside_x += data.vx[i] * data.m[i];
                momentum_inside_y += data.vy[i] * data.m[i];
                momentum_inside_z += data.vz[i] * data.m[i];
                mass_inside += data.m[i];
                n_inside++;
            }
        }
        momentum_inside_x *= 2.;
        momentum_inside_y *= 2.;
        momentum_inside_z *= 2.;
        mass_inside *= 2.;
        n_inside *= 2;
    }

    size_t countInside(size_t start, size_t end)
    {
        size_t count = 0;
        for (size_t i = start; i < end; i++)
        {
            T x    = data.x[i];
            T y    = data.y[i];
            T z    = data.z[i];
            T dist = std::sqrt(x * x + y * y + z * z);
            if (dist < inner_limit) count++;
        }
        return count;
    }

    size_t countInsideLocal() { return countInside(domain_ptr->startIndex(), domain_ptr->endIndex()); }

    size_t countInsideWithHalos() { return countInside(0, data.x.size()); }

    void sync()
    {
        std::vector<double>   s1, s2;
        std::vector<uint64_t> s3;
        std::vector<float>    s7, s8;

        domain_ptr->sync(data.keys, data.x, data.y, data.z, data.h, std::tie(data.m, data.vx, data.vy, data.vz),
                         std::tie(s1, s2, s3, s7, s8));
    }

    //! @brief Test if the particles that are marked for removal are inside the accretion radius
    void testComputeAccretionCondition()
    {
        for (size_t i = domain_ptr->startIndex(); i < domain_ptr->endIndex(); i++)
        {
            T x    = data.x[i];
            T y    = data.y[i];
            T z    = data.z[i];
            T dist = std::sqrt(x * x + y * y + z * z);
            if (data.keys[i] == 0) { EXPECT_TRUE(dist >= inner_limit); }
            else if (data.keys[i] == 1) { EXPECT_TRUE(dist < inner_limit); }
        }
    }

    void testMoveAccretedToEnd()
    {
        for (size_t i = domain_ptr->startIndex(); i < domain_ptr->endIndex() - star.n_accreted; i++)
        {
            T x    = data.x[i];
            T y    = data.y[i];
            T z    = data.z[i];
            T dist = std::sqrt(x * x + y * y + z * z);
            EXPECT_TRUE(dist >= inner_limit);
        }

        for (size_t i = domain_ptr->endIndex() - star.n_accreted; i < domain_ptr->endIndex(); i++)
        {
            T x    = data.x[i];
            T y    = data.y[i];
            T z    = data.z[i];
            T dist = std::sqrt(x * x + y * y + z * z);
            EXPECT_TRUE(dist < inner_limit);
        }
    }

    void testSumAccretedMassAndMomentum()
    {
        T momentum_sum_x = 0.;
        T momentum_sum_y = 0.;
        T momentum_sum_z = 0.;
        T m              = 0.;
        for (size_t i = domain_ptr->endIndex() - star.n_accreted; i < domain_ptr->endIndex(); i++)
        {

            momentum_sum_x += data.vx[i] * data.m[i];
            momentum_sum_y += data.vy[i] * data.m[i];
            momentum_sum_z += data.vz[i] * data.m[i];
            m += data.m[i];
        }
        EXPECT_NEAR(m, star.m_accreted_local, 1e-5);
        EXPECT_NEAR(momentum_sum_x, star.p_accreted_local[0], 1e-4);
        EXPECT_NEAR(momentum_sum_y, star.p_accreted_local[1], 1e-4);
        EXPECT_NEAR(momentum_sum_z, star.p_accreted_local[2], 1e-4);
    }
    void testExchangeAndAccreteOnStar()
    {
        // Check if mass of star = mass of accreted particles
        EXPECT_NEAR(star.m, 1. + mass_inside, 1e-9);
        // Check momentum of star
        EXPECT_NEAR(star.position_m1[0] / dt * star.m, momentum_inside_x, 1e-4);
        EXPECT_NEAR(star.position_m1[1] / dt * star.m, momentum_inside_y, 1e-4);
        EXPECT_NEAR(star.position_m1[2] / dt * star.m, momentum_inside_z, 1e-4);
    }
};

TEST_F(AccretionTest, testAccretion)
{
    if (numRanks != 2) throw std::runtime_error("Must be excuted with two ranks");

    initData();

    sync();

    const size_t first = domain_ptr->startIndex();
    const size_t last  = domain_ptr->endIndex();

    std::fill(data.keys.begin() + first, data.keys.end() + last, 0);

    planet::computeAccretionCondition(first, last, data, star);
    testComputeAccretionCondition();

    planet::computeNewOrder<ConservedFields, DependentFields>(first, last, data, star);
    planet::applyNewOrder<ConservedFields, DependentFields>(first, last, data, star);
    testMoveAccretedToEnd();

    planet::sumAccretedMassAndMomentum(first, last, data, star);
    testSumAccretedMassAndMomentum();

    planet::exchangeAndAccreteOnStar(star, dt, rank);
    testExchangeAndAccreteOnStar();

    domain_ptr->setEndIndex(last - star.n_accreted);

    sync();

    size_t n_inside            = countInsideLocal();
    size_t n_inside_with_halos = countInsideWithHalos();

    EXPECT_TRUE(n_inside == 0);
    EXPECT_TRUE(n_inside_with_halos == 0);
}