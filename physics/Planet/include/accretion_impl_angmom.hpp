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

template<typename Dataset, typename StarData>
// TODO: adjust function to remove particle based on threshold angular momentum 
// - calculate relative angular momentum of particle-star 
// - compare to threshold value 
void computeAccretionConditionImplAngMom(size_t first, size_t last, Dataset& d, StarData& star)
{
    const double star_size2 = star.inner_size * star.inner_size;
    const double R_Bondi = (2 * d.g * star.m) / (d.c * d.c); // bondi radius if want to use for sth
    const double J_circ2 = 2 * d.g * star.m *  star.inner_size; // angular momentum limit based on circular orbit 
    const double frac = 0.8; // fraction of escape angular momentum to consider
    const double J_esc2 = frac * J_circ2; // based on escape velocity

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
        // calculate specific angular momentum of particle i relative to star (cross product r and v)
        const double Jx = (dy*d.vz[i] - dz*d.vy[i]);
        const double Jy = (dz*d.vx[i] - dx*d.vz[i]);
        const double Jz = (dx*d.vy[i] - dy*d.vx[i]);
        const double J2 = Jx * Jx + Jy * Jy + Jz * Jz;

        // check if angular momentum is below limit 
        if (J2 < J_circ2) { remove_and_sum(i, accr_mass, accr_mom, n_accreted); }
        else if (d.h[i] > star.removal_limit_h) { remove_and_sum(i, removed_mass, removed_mom, n_removed); }
    }

    star.m_accreted_local    = accr_mass;
    star.p_accreted_local[0] = accr_mom[0];
    star.p_accreted_local[1] = accr_mom[1];
    star.p_accreted_local[2] = accr_mom[2];
    star.n_removed_local     = n_removed;

    star.m_removed_local    = removed_mass;
    star.p_removed_local[0] = removed_mom[0];
    star.p_removed_local[1] = removed_mom[1];
    star.p_removed_local[2] = removed_mom[2];
    star.n_accreted_local   = n_accreted;
}

} // namespace planet
