//
// Created by Noah Kubli on 12.03.2024.
//
#include <cub/cub.cuh>
#include <thrust/device_vector.h>
#include <thrust/partition.h>
#include <thrust/sequence.h>

#include "cstone/cuda/cuda_utils.cuh"
#include "cstone/findneighbors.hpp"
#include "cstone/traversal/find_neighbors.cuh"
#include "sph/util/device_math.cuh"

#include "cstone/sfc/box.hpp"
#include "cstone/tree/definitions.h"

#include "accretion_gpu.hpp"
#include "cuda_runtime.h"

static __device__ double dev_accr_mass;
static __device__ double dev_accr_mom_x;
static __device__ double dev_accr_mom_y;
static __device__ double dev_accr_mom_z;

template<typename T1, typename Th, typename Tremove, typename T2, typename Tm, typename Tv>
__global__ void computeAccretionConditionKernel(size_t first, size_t last, const T1* x, const T1* y, const T1* z,
                                                const Th* h, Tremove* remove, const Tm* m, const Tv* vx, const Tv* vy,
                                                const Tv* vz, T2 star_x, T2 star_y, T2 star_z, T2 star_size2,
                                                T2 removal_limit_h)
{
    cstone::LocalIndex i = first + blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= last) {}
    else
    {
        const double dx    = x[i] - star_x;
        const double dy    = y[i] - star_y;
        const double dz    = z[i] - star_z;
        const double dist2 = dx * dx + dy * dy + dz * dz;

        double accreted_mass{};
        double accreted_momentum_x;
        double accreted_momentum_y;
        double accreted_momentum_z;

        if (dist2 < star_size2)
        {
            remove[i]           = cstone::removeKey<Tremove>::value;
            accreted_mass       = m[i];
            accreted_momentum_x = m[i] * vx[i];
            accreted_momentum_y = m[i] * vy[i];
            accreted_momentum_z = m[i] * vz[i];

            /*remove[i] = 1;*/
        } // Accrete on star
        else if (h[i] > removal_limit_h)
        {
            remove[i] = cstone::removeKey<Tremove>::value; /*remove[i] = 2;*/
        }                                                  // Remove from system

        typedef cub::BlockReduce<Tm, TravConfig::numThreads> BlockReduceTm;
        __shared__ typename BlockReduceTm::TempStorage         temp_storage_m;
        BlockReduceTm                                          reduce_tm(temp_storage_m);
        Tm                                                   bs_accr_m = reduce_tm.Reduce(accreted_mass, cub::Sum());

        typedef cub::BlockReduce<Tv, TravConfig::numThreads> BlockReduceTv;
        __shared__ typename BlockReduceTv::TempStorage         temp_storage_px;
        __shared__ typename BlockReduceTv::TempStorage         temp_storage_py;
        __shared__ typename BlockReduceTv::TempStorage         temp_storage_pz;
        BlockReduceTv reduce_tvx(temp_storage_px);
        BlockReduceTv reduce_tvy(temp_storage_py);
        BlockReduceTv reduce_tvz(temp_storage_pz);

        Tv bs_accr_px = reduce_tvx.Reduce(accreted_momentum_x, cub::Sum());
        Tv bs_accr_py = reduce_tvy.Reduce(accreted_momentum_y, cub::Sum());
        Tv bs_accr_pz = reduce_tvz.Reduce(accreted_momentum_z, cub::Sum());

        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(&dev_accr_mass, bs_accr_m);
            atomicAdd(&dev_accr_mom_x, bs_accr_px);
            atomicAdd(&dev_accr_mom_y, bs_accr_py);
            atomicAdd(&dev_accr_mom_z, bs_accr_pz);
        }
    }
}

struct debug_zero
{
    __device__ bool operator()(size_t x) const { return x == 1; }
};

template<typename T1, typename Th, typename Tremove, typename T2, typename Tm, typename Tv>
void computeAccretionConditionGPU(size_t first, size_t last, const T1* x, const T1* y, const T1* z, const Th* h,
                                  Tremove* remove, const Tm* m, const Tv* vx, const Tv* vy, const Tv* vz,
                                  const T2* spos, T2 star_size, T2 removal_limit_h, Tm& m_accr, Tv& vx_accr,
                                  Tv& vy_accr, Tv& vz_accr)
{
    cstone::LocalIndex numParticles = last - first;
    unsigned           numThreads   = 256;
    unsigned           numBlocks    = (numParticles + numThreads - 1) / numThreads;

    static __device__ double dev_accr_mass;
    static __device__ double dev_accr_mom_x;
    static __device__ double dev_accr_mom_y;
    static __device__ double dev_accr_mom_z;

    double zero = 0.;
    cudaMemcpyToSymbol(dev_accr_mass, &zero, sizeof(zero));
    cudaMemcpyToSymbol(dev_accr_mom_x, &zero, sizeof(zero));
    cudaMemcpyToSymbol(dev_accr_mom_y, &zero, sizeof(zero));
    cudaMemcpyToSymbol(dev_accr_mom_z, &zero, sizeof(zero));
    computeAccretionConditionKernel<<<numBlocks, numThreads>>>(first, last, x, y, z, h, remove, m, vx, vy, vz, spos[0],
                                                               spos[1], spos[2], star_size * star_size,
                                                               removal_limit_h);
    checkGpuErrors(cudaGetLastError());
    checkGpuErrors(cudaDeviceSynchronize());

    double m_accr_ret;
    double px_accr_ret;
    double py_accr_ret;
    double pz_accr_ret;

    cudaMemcpyFromSymbol(&m_accr_ret, dev_accr_mass, sizeof(m_accr_ret));
    cudaMemcpyFromSymbol(&px_accr_ret, dev_accr_mom_x, sizeof(px_accr_ret));
    cudaMemcpyFromSymbol(&py_accr_ret, dev_accr_mom_y, sizeof(py_accr_ret));
    cudaMemcpyFromSymbol(&pz_accr_ret, dev_accr_mom_z, sizeof(pz_accr_ret));
    m_accr  = m_accr_ret;
    vx_accr = px_accr_ret;
    vy_accr = py_accr_ret;
    vz_accr = pz_accr_ret;

    // size_t nrem = thrust::count_if(thrust::device, remove + first, remove + last, debug_zero{});
    // printf("computeAccretionConditionGPU remove : %u\n", nrem);
}

template void computeAccretionConditionGPU(size_t, size_t, const double*, const double*, const double*, const float*,
                                           uint64_t*, const double*, const double*, const double*, const double*,
                                           const double*, double, double, double&, double&, double&, double&);

template void computeAccretionConditionGPU(size_t, size_t, const double*, const double*, const double*, const double*,
                                           uint64_t*, const double*, const double*, const double*, const double*,
                                           const double*, double, double, double&, double&, double&, double&);
template<typename T>
struct KeepParticle
{
    const T*        arr;
    __device__ bool operator()(const size_t& k) { return (arr[k] == 0); }
};

template<typename T>
struct AccreteParticle
{
    const T*        arr;
    __device__ bool operator()(const size_t& k) { return (arr[k] == 1); }
};

template<typename Tremove>
void computeNewOrderGPU(size_t first, size_t last, Tremove* remove, size_t* n_accreted, size_t* n_removed)
{
    thrust::device_vector<size_t> index(last - first);
    thrust::sequence(index.begin(), index.end(), first);

    const auto begin_accreted = thrust::stable_partition(index.begin(), index.end(), KeepParticle<Tremove>{remove});

    const auto begin_removed = thrust::stable_partition(begin_accreted, index.end(), AccreteParticle<Tremove>{remove});
    *n_accreted              = thrust::distance(begin_accreted, begin_removed);
    *n_removed               = thrust::distance(begin_removed, index.end());

    thrust::copy(thrust::device, index.begin(), index.end(), remove + first);
    checkGpuErrors(cudaDeviceSynchronize());
}

template void computeNewOrderGPU(size_t, size_t, size_t*, size_t*, size_t*);

template<typename T>
struct index_access
{
    const T*     x;
    __device__ T operator()(const size_t& k) { return (x[k]); }
};

template<typename T, typename Torder>
void applyNewOrderGPU(size_t first, size_t last, T* x, T* scratch, Torder* order)
{
    thrust::transform(thrust::device, order + first, order + last, scratch + first, index_access<T>{x});
    thrust::copy(thrust::device, scratch + first, scratch + last, x + first);
    checkGpuErrors(cudaDeviceSynchronize());
}

template void applyNewOrderGPU(size_t, size_t, double*, double*, size_t*);
template void applyNewOrderGPU(size_t, size_t, float*, float*, size_t*);
template void applyNewOrderGPU(size_t, size_t, unsigned long*, unsigned long*, size_t*);
template void applyNewOrderGPU(size_t, size_t, unsigned int*, unsigned int*, size_t*);

template<typename Tv, typename Tm, typename Tstar>
void sumMassAndMomentumGPU(size_t first, size_t last, const Tv* vx, const Tv* vy, const Tv* vz, const Tm* m,
                           Tv* scratch, Tstar* m_sum, Tstar* p_sum)
{
    if (first == last)
    {
        *m_sum   = 0.;
        p_sum[0] = 0.;
        p_sum[1] = 0.;
        p_sum[2] = 0.;
        return;
    }

    thrust::transform(thrust::device, vx + first, vx + last, m + first, scratch + first, thrust::multiplies<Tv>{});
    p_sum[0] = thrust::reduce(thrust::device, scratch + first, scratch + last, Tv{0.}, thrust::plus<Tv>{});

    thrust::transform(thrust::device, vy + first, vy + last, m + first, scratch + first, thrust::multiplies<Tv>{});
    p_sum[1] = thrust::reduce(thrust::device, scratch + first, scratch + last, Tv{0.}, thrust::plus<Tv>{});

    thrust::transform(thrust::device, vz + first, vz + last, m + first, scratch + first, thrust::multiplies<Tv>{});
    p_sum[2] = thrust::reduce(thrust::device, scratch + first, scratch + last, Tv{0.}, thrust::plus<Tv>{});

    *m_sum = thrust::reduce(thrust::device, m + first, m + last, Tm{0.}, thrust::plus<Tm>{});
}

template void sumMassAndMomentumGPU(size_t, size_t, const float*, const float*, const float*, const float*, float*,
                                    double*, double*);
template void sumMassAndMomentumGPU(size_t, size_t, const double*, const double*, const double*, const double*, double*,
                                    double*, double*);