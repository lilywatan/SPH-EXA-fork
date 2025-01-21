//
// Created by Noah Kubli on 15.03.2024.
//

#pragma once

#include <cassert>
#include <mpi.h>
#include "cstone/primitives/mpi_wrappers.hpp"
#include "cstone/fields/field_get.hpp"
#include "cstone/tree/accel_switch.hpp"
#include "cstone/cuda/cuda_stubs.h"
#include "cstone/util/type_list.hpp"

#include "sph/particles_data.hpp"

#include "accretion_impl_angmom.hpp"
#include "accretion_gpu.hpp"
#include "fieldListExclude.hpp"

namespace planet
{

//! @brief Flag particles for removal. Overwrites keys.
template<typename Dataset, typename StarData>
void computeAccretionConditionAngMom(size_t first, size_t last, Dataset& d, StarData& star)
{
    if constexpr (cstone::HaveGpu<typename Dataset::AcceleratorType>{})
    {
        computeAccretionConditionGPU(first, last, d, star);
    }
    else { computeAccretionConditionImplAngMom(first, last, d, star); }
}

} // namespace planet
