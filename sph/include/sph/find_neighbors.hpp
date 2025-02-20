#pragma once

#include "cstone/findneighbors.hpp"

namespace sph
{

using cstone::LocalIndex;


template<class Tc, class T, class KeyType>
void findNeighborsSph(const Tc* x, const Tc* y, const Tc* z, T* h, LocalIndex firstId, LocalIndex lastId,
                      const cstone::Box<Tc>& box, const cstone::OctreeNsView<Tc, KeyType>& treeView, unsigned ng0,
                      unsigned ngmax, unsigned ngmin, LocalIndex* neighbors, unsigned* nc)
{
    LocalIndex numWork = lastId - firstId;

    size_t        numFails     = 0;
    constexpr int maxIteration = 110;

#pragma omp parallel for reduction(+ : numFails)
    for (LocalIndex i = 0; i < numWork; ++i)
    {
        LocalIndex id    = i + firstId;
        unsigned   ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);

        T   h_upper(box.maxExtent());
        T   h_lower{0.};
        int iteration = 0;
        while ((ngmin > ncSph || (ncSph - 1) > ngmax) && iteration++ < maxIteration)
        {
            h_upper = (ncSph - 1) > ngmax ? h[i] : h_upper;
            h_lower = ngmin > ncSph ? h[i] : h_lower;
            if (iteration < 10)
            {
                // Dampen updateH by weighting with proposed smoothing lengths of past iterations
                h[id] = (updateH(ng0, ncSph, h[id]) + h[id] * iteration) / static_cast<T>(iteration + 1);
            }
            else
            {
                // Bisect algorithm
                h[i] = (h_upper + h_lower) / 2.;
            }
            ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);
        }

        numFails += (iteration >= maxIteration);

        nc[i] = ncSph;
    }

    if (numFails)
    {
        std::cout << "Coupled h-neighbor count updated failed to converge for " << numFails << " particles"
                  << std::endl;
    }
}


template<class Tc, class T, class KeyType>
void findNeighborsSph_debug(const Tc* x, const Tc* y, const Tc* z, T* h, LocalIndex firstId, LocalIndex lastId,
                      const cstone::Box<Tc>& box, const cstone::OctreeNsView<Tc, KeyType>& treeView, unsigned ng0,
                      unsigned ngmax, unsigned ngmin, LocalIndex* neighbors, unsigned* nc)
{
    LocalIndex numWork = lastId - firstId;
    size_t     numFails = 0;
    constexpr int maxIteration = 110;

    // Open summary debug file
    std::ofstream debugFile("neighbor_summary.log", std::ios::out);
    if (!debugFile) {
        std::cerr << "Error opening debug file!\n";
        return;
    }

    debugFile << "Iteration, Min_h, Max_h, Avg_h, Min_nc, Max_nc, Avg_nc, Convergence_Failures\n";

#pragma omp parallel for reduction(+ : numFails)
    for (LocalIndex i = 0; i < numWork; ++i)
    {
        LocalIndex id = i + firstId;
        unsigned   ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);

        T   h_upper(box.maxExtent());
        T   h_lower{0.};
        int iteration = 0;

        // Initialize statistics
        T min_h = std::numeric_limits<T>::max();
        T max_h = std::numeric_limits<T>::lowest();
        T sum_h = 0;
        unsigned min_nc = std::numeric_limits<unsigned>::max();
        unsigned max_nc = 0;
        unsigned sum_nc = 0;
        int failedParticles = 0;

        while ((ngmin > ncSph || (ncSph - 1) > ngmax) && iteration++ < maxIteration)
        {
            h_upper = (ncSph - 1) > ngmax ? h[id] : h_upper;
            h_lower = ngmin > ncSph ? h[id] : h_lower;

            if (iteration < 10)
            {
                h[id] = (updateH(ng0, ncSph, h[id]) + h[id] * iteration) / static_cast<T>(iteration + 1);
            }
            else
            {
                h[id] = (h_upper + h_lower) / 2.;
            }

            ncSph = 1 + findNeighbors(id, x, y, z, h, treeView, box, ngmax, neighbors + i * ngmax);

            // Update statistics
            min_h = std::min(min_h, h[id]);
            max_h = std::max(max_h, h[id]);
            sum_h += h[id];

            min_nc = std::min(min_nc, ncSph);
            max_nc = std::max(max_nc, ncSph);
            sum_nc += ncSph;
        }

        numFails += (iteration >= maxIteration);
        if (iteration >= maxIteration) {
            failedParticles++;
        }
        nc[i] = ncSph;

        // Safely compute averages, avoiding division by zero
#pragma omp critical
        {
            T avg_h = (iteration > 0) ? (sum_h / iteration) : h[id];
            unsigned avg_nc = (iteration > 0) ? (sum_nc / iteration) : ncSph;
            debugFile << iteration << ", " 
                      << min_h << ", " << max_h << ", " << avg_h << ", "
                      << min_nc << ", " << max_nc << ", " << avg_nc << ", "
                      << failedParticles << "\n";
        }
    }

    if (numFails) {
        debugFile << "WARNING: " << numFails << " particles failed to converge!\n";
        std::cout << "Coupled h-neighbor count updated failed to converge for " << numFails << " particles"
                  << std::endl;
    }

    debugFile.close();
}




//! @brief perform neighbor search together with updating the smoothing lengths
template<class T, class Dataset>
void findNeighborsSfc(size_t startIndex, size_t endIndex, Dataset& d, const cstone::Box<T>& box)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{}) { return; }

    if (d.ng0 > d.ngmax) { throw std::runtime_error("ng0 should be smaller than ngmax\n"); }
    if (d.ngmin > d.ng0) { throw std::runtime_error("ngmin should be smaller than ng0\n"); }


    findNeighborsSph_debug(d.x.data(), d.y.data(), d.z.data(), d.h.data(), startIndex, endIndex, box, d.treeView, d.ng0,
                     d.ngmax, d.ngmin, d.neighbors.data(), d.nc.data() + startIndex);
}

} // namespace sph
