#pragma once

#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>
#include <iostream>
#include "cstone/tree/definitions.h"

namespace planet
{
// function to take combine all locally accreted particles (at boundary of disk) 
// and compute the new masses of the disk & star based on the theoretical mass accretion rate
template<typename Dataset, typename DiskData, typename StarData>
void SubGridDiskAccreteOnStar(Dataset& d, DiskData& disk, StarData& star, double dt, int rank)
{
    // adjust calculations for sigma & mass accretion rate
    double m_accreted_global{};
    double h_accreted_global{};
    double r_disk_global{};
    double c_boundary_global_avg{};
    double n_accreted_global{};
    double r0_global{};
    double sigma0_global{};
    double rho_boundary_global_avg{};
    double T_boundary_global_avg{};
    double n_boundary_global{};
    double H2_boundary_global{};

    // Reductions on disk parameters
    MPI_Reduce(&disk.m_accreted_local, &m_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.h_accreted_local, &h_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.r_accreted_local, &r_disk_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.n_accreted_local, &n_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.c_boundary_local, &c_boundary_global_avg, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.r0_local, &r0_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.rho_boundary_local, &rho_boundary_global_avg, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.T_boundary_local, &T_boundary_global_avg, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.sigma0_local, &sigma0_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.n_boundary_local, &n_boundary_global, 1, MpiType<size_t>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.H2_boundary_local, &H2_boundary_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);

    // function to calculate mass accretion rate
    // nu is set disk kinematic viscosity 
    // r is the radius at which to evaluate the function 
    auto M_dot = [&star, &disk](double r)
    {
        double nu_0 = disk.alpha * disk.c_boundary * disk.Hr_boundary * disk.r0;
        double nu = disk.alpha * disk.c_boundary * disk.Hr_boundary * r;
        double star_r = 1 - std::sqrt(star.inner_size / r);
        double star_r0 = 1 - std::sqrt(star.inner_size / disk.r0);
        double sigma = disk.sigma0 * ((nu_0 * star_r) / (nu * star_r0));
        return (3 * M_PI * nu * disk.sigma0 * sigma) / star_r;
    };

    if (rank == 0){
        // change the disk.c thing
        // calculate rms of parameters
        h_accreted_global = std::sqrt(h_accreted_global / n_accreted_global);
        r_disk_global = std::sqrt(r_disk_global / n_accreted_global);
        c_boundary_global_avg = std::sqrt(c_boundary_global_avg / n_boundary_global);
        rho_boundary_global_avg = std::sqrt(rho_boundary_global_avg / n_boundary_global);
        T_boundary_global_avg = std::sqrt(T_boundary_global_avg / n_boundary_global);
        sigma0_global = std::sqrt(sigma0_global / n_boundary_global);
        r0_global = std::sqrt(r0_global / n_boundary_global);

        // set disk parameters to new values
        disk.h = h_accreted_global;
        disk.r = r_disk_global;
        disk.c_boundary = c_boundary_global_avg;
        disk.rho_boundary = rho_boundary_global_avg;
        disk.T_boundary = T_boundary_global_avg;
        disk.sigma0 = sigma0_global;
        disk.r0 = r0_global;
        disk.Hr_boundary =  std::sqrt(H2_boundary_global / m_accreted_global)/r0_global;

        // calculate mass accretion rate and new masses of star and disk 
        double m_star_new = star.m + M_dot(disk.r) * dt; 
        double m_disk_new = disk.m + m_accreted_global - M_dot(disk.r) * dt; 
        star.m = m_star_new;
        disk.m = m_disk_new;
    }

    MPI_Bcast(&star.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.r, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.h, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.c_boundary, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.rho_boundary, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.T_boundary, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.sigma0, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.r0, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.Hr_boundary, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);


}

} // namespace planet