//
// compute the subgrid disk structure 
// Created by  on 11.01.2024
//

// TO-DO : implement the surface density calculation at the boundary of the disk 
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
void SubGridDiskBoundary(size_t first, size_t last, Dataset& d, DiskData& disk, StarData& star)
{

    double boundary_density{};
    double boundary_temperature{}; 
    double boundary_sigma{};
    double boundary_mass{};
    size_t n_boundary_particles{};

    // sum density, temperature & mass of particles at the boundary of the disk, count number of particles
    // squared & multiplied by mass for mass weighted mean calculation
    auto add_parameters = [&d](size_t i, double& rho_sum, double& temp_sum, double& sigma_sum, double& mass_sum, size_t& n_boundary) 
    {
        d.keys[i] = cstone::removeKey<typename Dataset::KeyType>::value;
        rho_sum += d.rho[i]*d.rho[i]*d.m[i];
        temp_sum += d.temp[i]*d.temp[i]*d.m[i];
        mass_sum += d.m[i];
        n_boundary++;
    };

    // kernel function for surface density estimation
    auto W = [](double r, double h) {
        double q = r/h; 
        double sigma = 10 / (7 * M_PI); 
        if (0 <= q && q <= 1){
            return sigma / (h*h) * (1 - 1.5*q*q + 0.75*q*q*q);
        }
        else if (1 < q && q <= 2){
            return sigma / (h*h) * 0.25 * (2 - q)*(2 - q)*(2 - q);
        }
        else {
            return 0;
        }
    };

#pragma omp parallel for reduction(+ : boundary_density, boundary_temperature, boundary_sigma, boundary_mass, n_boundary_particles)
    for (size_t i = first; i < last; i++)
    {
        const double dx    = d.x[i] - star.position[0];
        const double dy    = d.y[i] - star.position[1];
        const double dz    = d.z[i] - star.position[2];
        const double dist2 = dx * dx + dy * dy + dz * dz;

        // TO-DO: implement some sort of radial boundary (interval) -> 2-3 times h currently 
        if (dist2 > 2 * d.h[i] && dist2 < 3 * d.h[i]) {
            add_parameters(i, boundary_density, boundary_temperature, boundary_sigma, boundary_mass, n_boundary_particles); 
            
        }
    }

    // mass weighted mean of density and temperature at disk boundary 
    if (boundary_mass > 0) {
        disk.rho_boundary_local = std::sqrt(boundary_density / boundary_mass);
        disk.temp_boundary_local = std::sqrt(boundary_temperature / boundary_mass);
    } 
    // if no mass then set to 0 
    else {
        disk.rho_boundary_local = 0;
        disk.temp_boundary_local = 0;
    }

    disk.n_boundary_local = n_boundary_particles;

}

} // namespace planet
