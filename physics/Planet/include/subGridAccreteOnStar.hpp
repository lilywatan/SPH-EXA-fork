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
    // finish editing this section -> 
    // add global reductions for disk parameters
    // adjust calculations for sigma & mass accretion rate
    double alpha = 0.1; 
    double m_accreted_global{};
    double c_global_avg{}
    double n_accreted_global{};
    double r0_global{};
    double sigma0_global{};
    double rho_boundary_global_avg{};
    double T_boundary_global_avg{};
    double n_boundary_global{};

    MPI_Reduce(&disk.m_accreted_local_subdisk, m_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.n_accreted_local, n_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    // reductions on averaging parameters -> rms calculation 
    MPI_Reduce(&disk.c_accreted_local, c_global_avg, 1, MpiType<double>{}, MPI_AVG, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.c_accreted_local, c_global_avg, 1, MpiType<double>{}, MPI_AVG, 0, MPI_COMM_WORLD);

    auto M_dot = [&star, &disk](double nu, double sigma)
    {
        double denominator = 1 - std::sqrt(star.inner_size / disk.r0);
        return (3 * M_PI * nu * disk.sigma_0 * sigma) / denominator;
    };

    auto sigma_norm = [] 
    {
        
    };
    if (rank == 0){
        // change the disk.c thing
        double mw_c_boundary = c_global_avg*c_global_avg*m_accreted_global;
        double mw_c_disk = disk.m * disk.c * disk.c;
        disk.c = std::sqrt((mw_c_boundary + mw_c_disk) / (m_accreted_global + disk.m));
        double nu = alpha * disk.c * disk.H_r * disk.r; // will be set? 
        double m_star_new = star.m + M_dot(nu, sigma) * dt; 
        double m_disk_new = disk.m + m_accreted_global - M_dot(nu, sigma) * dt; 
        star.m = m_star_new;
        disk.m = m_disk_new;
    }

    MPI_Bcast(&star.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.c, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);


}

} // namespace planet