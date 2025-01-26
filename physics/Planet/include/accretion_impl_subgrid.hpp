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
    double accr_c{}; 
    size_t n_accreted{};

    double boundary_mass{}; 
    double boundary_r0{}; 
    double boundary_rho{};
    double boundary_T{};
    double boundary_sigma0{};
    double n_boundary{};

    auto remove_and_sum = [&d](size_t i, double& mass_sum, size_t& n_sum, double& c_sum)
    {
        d.keys[i] = cstone::removeKey<typename Dataset::KeyType>::value;
        mass_sum += d.m[i];
        c_sum += d.c[i];
        n_sum++;
    };

    auto add_to_boundary = [&d](size_t i, double& mass_boundary, size_t& n_boundary, double& r0_boundary, double& rho_boundary, double& T_boundary, double dist2, double& sigma0_boundary)
    {
        mass_boundary += d.m[i];
        r0_boundary += dist2;
        rho_boundary += d.rho[i]*d.rho[i];
        T_boundary += d.u[i]*d.u[i];
        sigma0_boundary += d.m[i]/(M_PI*d.h[i]*d.h[i]);
        n_boundary++;
    };

#pragma omp parallel for reduction(+ : accr_mass) reduction(+ : accr_c) reduction(+ : n_accreted)              \
    reduction(+ : boundary_mass) reduction(+ : n_boundary) reduction(+ : boundary_r0) reduction(+ : boundary_rho)    \
    reduction(+ : boundary_T) reduction(+ : boundary_sigma0)
    for (size_t i = first; i < last; i++)
    {
        const double dx    = d.x[i] - star.position[0];
        const double dy    = d.y[i] - star.position[1];
        const double dz    = d.z[i] - star.position[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;

        // radial criterion based on smoothing length -> accrete onto disk: 
        if (dist2 < 2*d.h[i]) { remove_and_sum(i, accr_mass, n_accreted, accr_c,); }
        // radial criterion for boundary -> 2h < r < 3h & minimum number of neighbors 
        else if (dist2 > 2*d.h[i] && dist2 < 3*d.h[i] && d.nc >= 150) { add_to_boundary(i, boundary_mass, n_boundary, boundary_r0, boundary_rho, boundary_T, dist2, boundary_sigma0); }
        
        // Q: does the disk also need a removal limit? 
        //else if (d.h[i] > star.removal_limit_h) { remove_and_sum(i, removed_mass, removed_mom, n_removed); }
    }


    disk.m_accreted_local    = accr_mass;
    disk.c_accreted_local    = accr_c;
    disk.r0_local            = boundary_r0;
    disk.n_accreted_local    = n_accreted;

    disk.n_boundary_local    = n_boundary;
    disk.rho_boundary_local  = boundary_rho;
    disk.T_boundary_local    = boundary_T;
    disk.sigma0_local        = boundary_sigma0;

    
}

} // namespace planet
