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
    double removed_h{}; 
    double removed_r{}; 

    double boundary_c{}; 
    double boundary_mass{}; 
    double boundary_r0{}; 
    double boundary_rho{};
    double boundary_T{};
    double boundary_sigma0{};
    size_t n_boundary{};
    double boundary_H2{};

    double smoothing_length{};

    auto remove_and_sum = [&d](size_t i, double& mass_sum, double(&mom_sum)[3], size_t& n_sum, double& h_sum, double& r_sum, double dist2)
    {
        d.keys[i] = cstone::removeKey<typename Dataset::KeyType>::value;
        mass_sum += d.m[i];
        mom_sum[0] += d.m[i] * d.vx[i];
        mom_sum[1] += d.m[i] * d.vy[i];
        mom_sum[2] += d.m[i] * d.vz[i];
        h_sum += d.h[i]*d.h[i];
        r_sum += dist2;
        n_sum++;
    };

    auto add_to_boundary = [&d](size_t i, double& mass_boundary, size_t& n_boundary, double& r0_boundary, 
        double& rho_boundary, double& T_boundary, double dist2, double& sigma0_boundary, double& c_boundary, double& H2_boundary, double dz)
    {
        mass_boundary += d.m[i];
        r0_boundary += dist2;
        rho_boundary += d.rho[i]*d.rho[i];
        T_boundary += d.u[i]*d.u[i];
        c_boundary += d.c[i]*d.c[i];
        double sigma0 = d.m[i]/(M_PI*d.h[i]*d.h[i]);
        sigma0_boundary += sigma0 * sigma0;
        n_boundary++;
        H2_boundary += d.m[i] * dz * dz;
    };

#pragma omp parallel for reduction(+ : accr_mass) reduction(+ : accr_mom[ : 3]) reduction(+ : boundary_c) reduction(+ : n_accreted) \
    reduction(+ : boundary_mass) reduction(+ : n_boundary) reduction(+ : boundary_r0) reduction(+ : boundary_rho)    \
    reduction(+ : boundary_T) reduction(+ : boundary_sigma0) reduction(+ : boundary_H2) reduction(+: removed_h) \
    reduction(+ :  removed_r) reduction(+ : smoothing_length) 
    for (size_t i = first; i < last; i++)
    {
        const double dx    = d.x[i] - star.position[0];
        const double dy    = d.y[i] - star.position[1];
        const double dz    = d.z[i] - star.position[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;

        smoothing_length += d.h[i];
        const double min_h = 2; 
        const double max_h = 3;
        // radial criterion based on smoothing length -> accrete onto disk: 
        if (dist2 < (min_h*d.h[i])*(min_h*d.h[i])) { 
            remove_and_sum(i, accr_mass, accr_mom, n_accreted, removed_h, removed_r, dist2); }
        // radial criterion for boundary -> 2h < r < 3h 
        else if (dist2 > (min_h*d.h[i])*(min_h*d.h[i]) && dist2 < (max_h*d.h[i])*(max_h*d.h[i])) { add_to_boundary(i, boundary_mass, n_boundary, 
            boundary_r0, boundary_rho, boundary_T, dist2, boundary_sigma0, boundary_c, boundary_H2, dz); }
        
        /* if (dist2 < 2 * 2 && d.h[i] < 2.0) { 
            remove_and_sum(i, accr_mass, n_accreted, removed_h, removed_r, dist2); } 

        else if (dist2 > 2 * 2 && dist2 < 3 * 3 && d.h[i] < 2.0) { add_to_boundary(i, boundary_mass, n_boundary, 
            boundary_r0, boundary_rho, boundary_T, dist2, boundary_sigma0, boundary_c, boundary_H2, dz); }
            */
        //else if (d.h[i] > star.removal_limit_h) { remove_and_sum(i, removed_mass, removed_mom, n_removed); }
    }
    std::cout << "Smoothing length: " << smoothing_length / d.numParticlesGlobal << std::endl;

    disk.m_accreted_local    = accr_mass;
    disk.r0_local            = boundary_r0;
    disk.n_accreted_local    = n_accreted;
    disk.h_accreted_local    = removed_h;
    disk.r_accreted_local    = removed_r;
    star.p_accreted_local[0] = accr_mom[0];
    star.p_accreted_local[1] = accr_mom[1];
    star.p_accreted_local[2] = accr_mom[2];

    disk.c_boundary_local    = boundary_c;
    disk.n_boundary_local    = n_boundary;
    disk.rho_boundary_local  = boundary_rho;
    disk.T_boundary_local    = boundary_T;
    disk.sigma0_local        = boundary_sigma0;
    disk.H2_boundary_local   = boundary_H2;

    
}

} // namespace planet
