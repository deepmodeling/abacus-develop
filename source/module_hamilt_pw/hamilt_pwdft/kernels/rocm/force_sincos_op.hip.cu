#include "module_hamilt_pw/hamilt_pwdft/kernels/force_sincos_op.h"

#include <complex>
#include <thrust/complex.h>
#include <hip/hip_runtime.h>
#include <base/macros/macros.h>

#define THREADS_PER_BLOCK 256

namespace hamilt {

// HIP kernels for sincos loops
template <typename FPTYPE>
__global__ void cal_force_loc_sincos_kernel(
    const int nat,
    const int npw,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const int* iat2it,
    const FPTYPE* vloc_factors,
    const thrust::complex<FPTYPE>* aux,
    const FPTYPE scale_factor,
    FPTYPE* force)
{
    const FPTYPE TWO_PI = 2.0 * M_PI;
    
    const int iat = blockIdx.y;
    const int ig_start = blockIdx.x * blockDim.x + threadIdx.x;
    const int ig_stride = gridDim.x * blockDim.x;
    
    if (iat >= nat) return;
    
    // Load atom information to registers
    const int it = iat2it[iat];
    const FPTYPE tau_x = tau[iat * 3 + 0];
    const FPTYPE tau_y = tau[iat * 3 + 1];
    const FPTYPE tau_z = tau[iat * 3 + 2];
    
    // Local accumulation variables
    FPTYPE local_force_x = 0.0;
    FPTYPE local_force_y = 0.0;
    FPTYPE local_force_z = 0.0;
    
    // Grid-stride loop over G-vectors
    for (int ig = ig_start; ig < npw; ig += ig_stride) {
        // Calculate phase: 2π * (G · τ)
        const FPTYPE phase = TWO_PI * (
            gcar[ig * 3 + 0] * tau_x + 
            gcar[ig * 3 + 1] * tau_y + 
            gcar[ig * 3 + 2] * tau_z
        );
        
        // Use HIP intrinsic for sincos
        FPTYPE sinp, cosp;
        __sincos(phase, &sinp, &cosp);
        
        // Calculate force factor
        const FPTYPE factor = vloc_factors[ig] * 
                              (cosp * aux[ig].imag() + sinp * aux[ig].real());
        
        // Accumulate force contributions
        local_force_x += gcar[ig * 3 + 0] * factor;
        local_force_y += gcar[ig * 3 + 1] * factor;
        local_force_z += gcar[ig * 3 + 2] * factor;
    }
    
    // Atomic add to global memory
    atomicAdd(&force[iat * 3 + 0], local_force_x * scale_factor);
    atomicAdd(&force[iat * 3 + 1], local_force_y * scale_factor);
    atomicAdd(&force[iat * 3 + 2], local_force_z * scale_factor);
}

template <typename FPTYPE>
__global__ void cal_force_ew_sincos_kernel(
    const int nat,
    const int npw,
    const int ig_gge0,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const int* iat2it,
    const FPTYPE* it_facts,
    const thrust::complex<FPTYPE>* aux,
    FPTYPE* force)
{
    const FPTYPE TWO_PI = 2.0 * M_PI;
    
    const int iat = blockIdx.y;
    const int ig_start = blockIdx.x * blockDim.x + threadIdx.x;
    const int ig_stride = gridDim.x * blockDim.x;
    
    if (iat >= nat) return;
    
    // Load atom information to registers
    const int it = iat2it[iat];
    const FPTYPE tau_x = tau[iat * 3 + 0];
    const FPTYPE tau_y = tau[iat * 3 + 1];
    const FPTYPE tau_z = tau[iat * 3 + 2];
    const FPTYPE it_fact = it_facts[iat];
    
    // Local accumulation variables
    FPTYPE local_force_x = 0.0;
    FPTYPE local_force_y = 0.0;
    FPTYPE local_force_z = 0.0;
    
    // Grid-stride loop over G-vectors
    for (int ig = ig_start; ig < npw; ig += ig_stride) {
        // Skip G=0 term
        if (ig == ig_gge0) continue;
        
        // Calculate phase: 2π * (G · τ)
        const FPTYPE arg = TWO_PI * (
            gcar[ig * 3 + 0] * tau_x + 
            gcar[ig * 3 + 1] * tau_y + 
            gcar[ig * 3 + 2] * tau_z
        );
        
        // Use HIP intrinsic for sincos
        FPTYPE sinp, cosp;
        __sincos(arg, &sinp, &cosp);
        
        // Calculate Ewald sum contribution
        const FPTYPE sumnb = -cosp * aux[ig].imag() + sinp * aux[ig].real();
        
        // Accumulate force contributions
        local_force_x += gcar[ig * 3 + 0] * sumnb;
        local_force_y += gcar[ig * 3 + 1] * sumnb;
        local_force_z += gcar[ig * 3 + 2] * sumnb;
    }
    
    // Atomic add to global memory
    atomicAdd(&force[iat * 3 + 0], local_force_x * it_fact);
    atomicAdd(&force[iat * 3 + 1], local_force_y * it_fact);
    atomicAdd(&force[iat * 3 + 2], local_force_z * it_fact);
}

// GPU operator implementations
template <typename FPTYPE>
void cal_force_loc_sincos_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* ctx,
    const int& nat,
    const int& npw,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const int* iat2it,
    const FPTYPE* vloc_factors,
    const std::complex<FPTYPE>* aux,
    const FPTYPE& scale_factor,
    FPTYPE* force)
{
    // Calculate optimal grid configuration
    const int threads_per_block = THREADS_PER_BLOCK;
    const int max_blocks_x = min(32, (npw + threads_per_block - 1) / threads_per_block);
    
    dim3 grid(max_blocks_x, nat);
    dim3 block(threads_per_block);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_loc_sincos_kernel<FPTYPE>), grid, block, 0, 0,
        nat, npw,
        gcar, tau, iat2it, vloc_factors,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(aux),
        scale_factor, force
    );
    
    hipCheckOnDebug();
}

template <typename FPTYPE>
void cal_force_ew_sincos_op<FPTYPE, base_device::DEVICE_GPU>::operator()(
    const base_device::DEVICE_GPU* ctx,
    const int& nat,
    const int& npw,
    const int& ig_gge0,
    const FPTYPE* gcar,
    const FPTYPE* tau,
    const int* iat2it,
    const FPTYPE* it_facts,
    const std::complex<FPTYPE>* aux,
    FPTYPE* force)
{
    // Calculate optimal grid configuration
    const int threads_per_block = THREADS_PER_BLOCK;
    const int max_blocks_x = min(32, (npw + threads_per_block - 1) / threads_per_block);
    
    dim3 grid(max_blocks_x, nat);
    dim3 block(threads_per_block);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cal_force_ew_sincos_kernel<FPTYPE>), grid, block, 0, 0,
        nat, npw, ig_gge0,
        gcar, tau, iat2it, it_facts,
        reinterpret_cast<const thrust::complex<FPTYPE>*>(aux),
        force
    );
    
    hipCheckOnDebug();
}

template struct cal_force_loc_sincos_op<float, base_device::DEVICE_GPU>;
template struct cal_force_loc_sincos_op<double, base_device::DEVICE_GPU>;

template struct cal_force_ew_sincos_op<float, base_device::DEVICE_GPU>;
template struct cal_force_ew_sincos_op<double, base_device::DEVICE_GPU>;

}  // namespace hamilt 