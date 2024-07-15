#pragma once

#include <vector>

#include "cstone/util/tuple.hpp"

#include "kernels.hpp"

namespace sph
{

//! @brief returns the heat capacity for given mean molecular weight
template<class T1, class T2>
HOST_DEVICE_FUN constexpr T1 idealGasCv(T1 mui, T2 gamma)
{
    constexpr T1 R = 8.317e7;
    return R / mui / (gamma - T1(1));
}

/*! @brief Reduced version of Ideal gas EOS for internal energy
 *
 * @param u     internal energy
 * @param rho   baryonic density
 * @param mui   mean molecular weight
 * @param gamma adiabatic index
 *
 * This EOS is used for simple cases where we don't need the temperature.
 * Returns pressure, speed of sound
 */
template<class T1, class T2, class T3>
HOST_DEVICE_FUN auto idealGasEOS_u(T1 u, T2 rho, T3 gamma)
{
    using Tc = std::common_type_t<T1, T2, T3>;

    Tc tmp = u * (gamma - Tc(1));
    Tc p   = rho * tmp;
    Tc c   = std::sqrt(gamma * tmp);

    return util::tuple<Tc, Tc>{p, c};
}

template<class T1, class T2, class T3>
HOST_DEVICE_FUN auto idealGasEOSTemp(T1 temp, T2 rho, T3 mui, T1 gamma)
{
    using Tc = std::common_type_t<T1, T2, T3>;
    return idealGasEOS_u(idealGasCv(mui, gamma) * temp, rho, gamma);
    /*Tc tmp = idealGasCv(mui, gamma) * temp * (gamma - Tc(1));
    Tc p   = rho * tmp;
    Tc c   = std::sqrt(tmp);

    return util::tuple<Tc, Tc>{p, c};*/
}

/*! @brief Polytropic EOS for a 1.4 M_sun and 12.8 km neutron star
 *
 * @param rho  baryonic density
 *
 * Kpol is hardcoded for these NS characteristics and is not valid for
 * other NS masses and radius
 * Returns pressure, and speed of sound
 */
template<class T>
HOST_DEVICE_FUN auto polytropicEOS_neutronstar(T rho)
{
    constexpr T Kpol     = 2.246341237993810232e-10;
    constexpr T gammapol = 3.e0;

    T p = Kpol * std::pow(rho, gammapol);
    T c = std::sqrt(gammapol * p / rho);

    return util::tuple<T, T>{p, c};
}

/*! @brief Polytropic EOS interface for SPH where rho is computed on-the-fly
 *
 * @tparam Dataset
 * @param startIndex  index of first locally owned particle
 * @param endIndex    index of last locally owned particle
 * @param d           the dataset with the particle buffers
 */
template<typename Dataset>
void computeEOS_Polytropic_neutronstar(size_t startIndex, size_t endIndex, Dataset& d)
{
    const auto* kx = d.kx.data();
    const auto* xm = d.xm.data();
    const auto* m  = d.m.data();

    auto* p = d.p.data();
    auto* c = d.c.data();

#pragma omp parallel for schedule(static)
    for (size_t i = startIndex; i < endIndex; ++i)
    {
        auto rho             = kx[i] * m[i] / xm[i];
        std::tie(p[i], c[i]) = polytropicEOS_neutronstar(rho);
    }
}

/*! @brief general polytropic equation of state
 *
 * @param rho SPH density
 */
template<typename T1, typename T2, typename T3, typename T4>
HOST_DEVICE_FUN auto polytropicEOS(T1 Kpoly, T2 exp_poly, T3 gamma, T4 rho)
{
    using Tc = std::common_type_t<T1, T2, T3, T4>;

    Tc p = Kpoly * std::pow(rho, exp_poly);
    Tc c = std::sqrt(gamma * p / rho);

    return util::tuple<Tc, Tc>{p, c};
}

} // namespace sph
