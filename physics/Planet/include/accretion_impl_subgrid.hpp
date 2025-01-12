//
// adjusted ConditionImplAng function to remove particles based on angular momentum
// Created by Noah Kubli on 14.03.2024.
//

#pragma once

#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>
#include <iostream>
#include "cstone/tree/definitions.h"

namespace planet
{

template<typename Dataset, typename DiskData, typename StarData>
// TODO: implement accretion based on radial condition for subdisk 
void computeAccretionConditionImplSubGridDisk(size_t first, size_t last, Dataset& d, DiskData& disk, StarData& star)
{

    double accr_mass{};
    double accr_mom[3]{};
    size_t n_accreted{};

    double removed_mass{};
    double removed_mom[3]{};
    size_t n_removed{};

    auto remove_and_sum = [&d](size_t i, double& mass_sum, double(&mom_sum)[3], size_t& n_sum)
    {
        d.keys[i] = cstone::removeKey<typename Dataset::KeyType>::value;
        mass_sum += d.m[i];
        mom_sum[0] += d.m[i] * d.vx[i];
        mom_sum[1] += d.m[i] * d.vy[i];
        mom_sum[2] += d.m[i] * d.vz[i];
        n_sum++;
    };

#pragma omp parallel for reduction(+ : accr_mass) reduction(+ : accr_mom[ : 3]) reduction(+ : n_accreted)              \
    reduction(+ : n_removed)
    for (size_t i = first; i < last; i++)
    {
        const double dx    = d.x[i] - star.position[0];
        const double dy    = d.y[i] - star.position[1];
        const double dz    = d.z[i] - star.position[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;

        // radial criterion: 
        if (dist2 < disk.r_out) { remove_and_sum(i, accr_mass, accr_mom, n_accreted); }
        // Q: does the disk also need a removal limit? 
        else if (d.h[i] > star.removal_limit_h) { remove_and_sum(i, removed_mass, removed_mom, n_removed); }
    }

    // add these to Disk Data? depends on disk-star interaction
    disk.m_accreted_local    = accr_mass;
    disk.p_accreted_local[0] = accr_mom[0];
    disk.p_accreted_local[1] = accr_mom[1];
    disk.p_accreted_local[2] = accr_mom[2];
    disk.n_removed_local     = n_removed;

    disk.m_removed_local    = removed_mass;
    disk.p_removed_local[0] = removed_mom[0];
    disk.p_removed_local[1] = removed_mom[1];
    disk.p_removed_local[2] = removed_mom[2];
    disk.n_accreted_local   = n_accreted;
}

} // namespace planet
