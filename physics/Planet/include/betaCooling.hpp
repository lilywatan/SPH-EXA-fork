//
// Created by Noah Kubli on 17.04.2024.
//

#pragma once
#include <cmath>
#include "cstone/tree/accel_switch.hpp"
#include "betaCooling_gpu.hpp"
#include "sph/particles_data.hpp"
namespace planet
{

template<typename Tpos, typename Ts, typename Tdu, typename Trho, typename Tu>
void betaCoolingImpl(size_t first, size_t last, const Tpos* x, const Tpos* y, const Tpos* z, Tdu* du, Tu* u,
                     Ts star_mass, const Ts* star_pos, Ts beta, Tpos g, Trho* rho,
                     Trho cooling_rho_limit = 1.683e-3)
{
    double cooling_floor = 9.3e-6; // approx. 1 K;
    //2.73 K: u = 2.5e-5;

    //Changed if condition (not yet started)
    size_t n_below_floor{};
    size_t n_nan{};
    for (size_t i = first; i < last; i++)
    {
        if (rho[i] < cooling_rho_limit && u[i] > cooling_floor) {

            const double dx    = x[i] - star_pos[0];
            const double dy    = y[i] - star_pos[1];
            const double dz    = z[i] - star_pos[2];
            const double dist2 = dx * dx + dy * dy + dz * dz;
            const double dist  = std::sqrt(dist2);
            const double omega = std::sqrt(g * star_mass / (dist2 * dist));
            du[i] += -u[i] * omega / beta;
        }

        if (u[i] < cooling_floor)
        {
            u[i] = cooling_floor;
            du[i] = std::max(0., du[i]);
            n_below_floor++;
        }
        if (std::isnan(du[i]))
            {
            n_nan++;
        }

    }
    printf("n_below_floor: %zu\n", n_below_floor);
    if (n_nan > 0) printf("Have nan particles: %zu\n", n_nan);
}

template<typename Dataset, typename StarData>
void betaCooling(Dataset& d, size_t startIndex, size_t endIndex, const StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        transferToHost(d, startIndex, endIndex, {"x", "y", "z", "du", "u", "rho"});
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.du.data(), d.u.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), star.cooling_rho_limit);
        transferToDevice(d, startIndex, endIndex, {"du", "u"});

        /*betaCoolingGPU(startIndex, endIndex, rawPtr(d.devData.x), rawPtr(d.devData.y), rawPtr(d.devData.z),
                       rawPtr(d.devData.du), rawPtr(d.devData.u), star.m, star.position.data(), star.beta, d.g,
                       rawPtr(d.devData.rho), star.cooling_rho_limit);*/
    }
    else
    {
        betaCoolingImpl(startIndex, endIndex, d.x.data(), d.y.data(), d.z.data(), d.du.data(), d.u.data(), star.m,
                        star.position.data(), star.beta, d.g, d.rho.data(), star.cooling_rho_limit);
    }
}
} // namespace planet
