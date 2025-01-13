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
    double alpha = 0.1; 
    double m_accreted_global{};
    double c_global_avg{}; 
    MPI_Reduce(&disk.m_accreted_local_subdisk, m_accreted_global, 1, MpiType<double>{}, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&disk.sound_speed_boundary_local, c_global_avg, 1, MpiType<double>{}, MPI_AVG, 0, MPI_COMM_WORLD);

    auto M_dot = [&d](double nu, double sigma)
    {
        return 3 * M_PI * nu * sigma;
    };

    if (rank == 0){
        double mw_c_boundary = c_global_avg*c_global_avg*m_accreted_global;
        double mw_c_disk = disk.m * disk.c * disk.c;
        disk.c = std::sqrt((mw_c_boundary + mw_c_disk) / (m_accreted_global + disk.m));
        double nu = alpha * disk.c * disk.H; 
        double m_star_new = star.m + M_dot(nu, sigma) * dt; 
        double m_disk_new = disk.m + m_accreted_global - M_dot(nu, sigma) * dt; 
        star.m = m_star_new;
        disk.m = m_disk_new;
    }

    MPI_Bcast(&star.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);
    MPI_Bcast(&disk.m, 1, MpiType<double>{}, 0, MPI_COMM_WORLD);


}

} // namespace planet