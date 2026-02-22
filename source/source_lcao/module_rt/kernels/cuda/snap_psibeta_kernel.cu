/**
 * @file snap_psibeta_kernel.cu
 * @brief CUDA kernel implementation for <psi|beta> overlap integrals
 *
 * This file implements the GPU-accelerated numerical integration for computing
 * overlap integrals between atomic orbitals (psi) and non-local projectors (beta).
 * The implementation uses:
 * - Lebedev-Laikov quadrature (110 points) for angular integration
 * - Gauss-Legendre quadrature (140 points) for radial integration
 * - Templated spherical harmonics with compile-time L for optimization
 * - Warp-level shuffle reduction for efficient parallel summation
 */

#include "snap_psibeta_kernel.cuh"
#include "source_base/constants.h"
#include "source_base/math_integral.h"

#include <cstdio>
#include <vector>

namespace module_rt
{
namespace gpu
{

//=============================================================================
// Constant Memory - Integration Grids
//=============================================================================

// Lebedev-Laikov angular quadrature grid (110 points)
__constant__ double d_lebedev_x[ANGULAR_GRID_NUM]; ///< x-direction cosines
__constant__ double d_lebedev_y[ANGULAR_GRID_NUM]; ///< y-direction cosines
__constant__ double d_lebedev_z[ANGULAR_GRID_NUM]; ///< z-direction cosines
__constant__ double d_lebedev_w[ANGULAR_GRID_NUM]; ///< Angular integration weights

// Gauss-Legendre radial quadrature grid (140 points)
__constant__ double d_gl_x[RADIAL_GRID_NUM]; ///< Quadrature abscissae on [-1, 1]
__constant__ double d_gl_w[RADIAL_GRID_NUM]; ///< Quadrature weights

//=============================================================================
// Warp-Level Reduction
//=============================================================================

/**
 * @brief Warp-level sum reduction using shuffle instructions
 *
 * Performs a parallel reduction within a warp (32 threads) using __shfl_down_sync.
 * After this function, lane 0 contains the sum of all input values in the warp.
 *
 * @param val Input value from each thread
 * @return Sum across all threads in the warp (valid only in lane 0)
 */
__device__ __forceinline__ double warp_reduce_sum(double val)
{
    for (int offset = 16; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

//=============================================================================
// Main Kernel Implementation
//=============================================================================

/**
 * @brief Atom-level batch kernel for <psi|beta> overlap integrals
 *
 * Integration is performed using restructured loops for efficiency:
 * - Outer loop: angular points (each thread handles different angles)
 * - Inner loop: radial points (each thread accumulates all radii)
 *
 * This structure exploits the fact that Y_lm for the projector (ylm0) only
 * depends on the angular direction, not the radial distance, saving
 * RADIAL_GRID_NUM redundant ylm0 computations per angular point.
 */
__global__ void snap_psibeta_atom_batch_kernel(double3 R0,
                                               double3 A,
                                               const NeighborOrbitalData* __restrict__ neighbor_orbitals,
                                               const ProjectorData* __restrict__ projectors,
                                               const double* __restrict__ psi_radial,
                                               const double* __restrict__ beta_radial,
                                               const int* __restrict__ proj_m0_offset,
                                               int total_neighbor_orbitals,
                                               int nproj,
                                               int natomwfc,
                                               int nlm_dim,
                                               cuDoubleComplex* __restrict__ nlm_out)
{
    // Thread/block indices
    const int norb_idx = blockIdx.x; // Which (neighbor, orbital) pair
    const int proj_idx = blockIdx.y; // Which projector
    const int tid = threadIdx.x;

    // Early exit for out-of-bounds blocks
    if (norb_idx >= total_neighbor_orbitals || proj_idx >= nproj)
    {
        return;
    }

    //-------------------------------------------------------------------------
    // Load input data
    //-------------------------------------------------------------------------

    const NeighborOrbitalData& norb = neighbor_orbitals[norb_idx];
    const ProjectorData& proj = projectors[proj_idx];

    const double3 R1 = norb.R1;
    const int L1 = norb.L1;
    const int m1 = norb.m1;
    const int L0 = proj.L0;
    const int m0_offset = proj_m0_offset[proj_idx];

    // Skip if angular momentum exceeds supported limit
    if (L1 > MAX_L || L0 > MAX_L)
    {
        return;
    }

    //-------------------------------------------------------------------------
    // Compute geometry
    //-------------------------------------------------------------------------

    // Note: dR (R1 - R0) is computed inline as dRx/dRy/dRz in the integration loop

    // Orbital cutoff
    const double r1_max = norb.psi_rcut;

    // Integration range from projector radial grid
    const double r_min = proj.r_min;
    const double r_max = proj.r_max;
    const double xl = 0.5 * (r_max - r_min);    // Half-range for Gauss-Legendre
    const double xmean = 0.5 * (r_max + r_min); // Midpoint

    // Phase factor exp(i * A · R0)
    const double AdotR0 = A.x * R0.x + A.y * R0.y + A.z * R0.z;
    const cuDoubleComplex exp_iAR0 = cu_exp_i(AdotR0);

    //-------------------------------------------------------------------------
    // Shared memory for warp reduction
    //-------------------------------------------------------------------------

    constexpr int NUM_WARPS = BLOCK_SIZE / 32; // 128 / 32 = 4 warps
    __shared__ double s_temp_re[NUM_WARPS];
    __shared__ double s_temp_im[NUM_WARPS];

    //-------------------------------------------------------------------------
    // Initialize accumulators (per-thread registers)
    //-------------------------------------------------------------------------

    const int num_m0 = 2 * L0 + 1;

    double result_re[MAX_M0_SIZE];
    double result_im[MAX_M0_SIZE];
    double result_r_re[3][MAX_M0_SIZE]; // For current operator: x, y, z components
    double result_r_im[3][MAX_M0_SIZE];

    for (int m0 = 0; m0 < num_m0; m0++)
    {
        result_re[m0] = 0.0;
        result_im[m0] = 0.0;
        for (int d = 0; d < 3; d++)
        {
            result_r_re[d][m0] = 0.0;
            result_r_im[d][m0] = 0.0;
        }
    }

    //-------------------------------------------------------------------------
    // Main integration loop
    // Outer: angular points (parallelized across threads)
    // Inner: radial points (accumulated per thread)
    //-------------------------------------------------------------------------

    for (int ian = tid; ian < ANGULAR_GRID_NUM; ian += BLOCK_SIZE)
    {
        // Load angular grid point
        const double leb_x = d_lebedev_x[ian];
        const double leb_y = d_lebedev_y[ian];
        const double leb_z = d_lebedev_z[ian];
        const double w_ang = d_lebedev_w[ian];

        // Precompute Y_lm for projector (independent of radial distance)
        double ylm0[MAX_YLM_SIZE];
        ModuleBase::sph_harm(L0, leb_x, leb_y, leb_z, ylm0);
        const int offset_L0 = L0 * L0;

        // Precompute A · direction (for phase factor)
        const double A_dot_leb = A.x * leb_x + A.y * leb_y + A.z * leb_z;

        // Vector from R1 to R0 (for computing distance to orbital center)
        const double dRx = R0.x - R1.x;
        const double dRy = R0.y - R1.y;
        const double dRz = R0.z - R1.z;

        // Radial integration
#pragma unroll 4
        for (int ir = 0; ir < RADIAL_GRID_NUM; ir++)
        {
            // Transform Gauss-Legendre point from [-1,1] to [r_min, r_max]
            const double r_val = xmean + xl * d_gl_x[ir];
            const double w_rad = xl * d_gl_w[ir];

            // Integration point position relative to R0
            const double rx = r_val * leb_x;
            const double ry = r_val * leb_y;
            const double rz = r_val * leb_z;

            // Vector from R1 to integration point
            const double tx = rx + dRx;
            const double ty = ry + dRy;
            const double tz = rz + dRz;
            const double tnorm = sqrt(tx * tx + ty * ty + tz * tz);

            // Check if within orbital cutoff
            if (tnorm <= r1_max)
            {
                // Compute Y_lm for orbital (depends on direction from R1)
                double ylm1[MAX_YLM_SIZE];
                if (tnorm > 1e-10)
                {
                    const double inv_tnorm = 1.0 / tnorm;
                    ModuleBase::sph_harm(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, ylm1);
                }
                else
                {
                    ModuleBase::sph_harm(L1, 0.0, 0.0, 1.0, ylm1);
                }

                // Interpolate orbital radial function
                const double psi_val
                    = interpolate_radial_gpu(psi_radial + norb.psi_offset, norb.psi_mesh, 1.0 / norb.psi_dk, tnorm);

                // Interpolate projector radial function
                const double beta_val
                    = interpolate_radial_gpu(beta_radial + proj.beta_offset, proj.beta_mesh, 1.0 / proj.beta_dk, r_val);

                // Phase factor exp(i * A · r)
                const double phase = r_val * A_dot_leb;
                const cuDoubleComplex exp_iAr = cu_exp_i(phase);

                // Orbital Y_lm value
                const double ylm_L1_val = ylm1[L1 * L1 + m1];

                // Combined integration factor: Y_L1m1 * psi * beta * r * dr * dOmega
                const double factor = ylm_L1_val * psi_val * beta_val * r_val * w_rad * w_ang;
                const cuDoubleComplex common_factor = cu_mul_real(exp_iAr, factor);

                // Accumulate for all m0 components of projector
#pragma unroll
                for (int m0 = 0; m0 < num_m0; m0++)
                {
                    const double ylm0_val = ylm0[offset_L0 + m0];

                    result_re[m0] += common_factor.x * ylm0_val;
                    result_im[m0] += common_factor.y * ylm0_val;

                    // Current operator contribution (if requested)
                    if (nlm_dim == 4)
                    {
                        const double r_op_x = rx + R0.x;
                        const double r_op_y = ry + R0.y;
                        const double r_op_z = rz + R0.z;

                        result_r_re[0][m0] += common_factor.x * ylm0_val * r_op_x;
                        result_r_im[0][m0] += common_factor.y * ylm0_val * r_op_x;
                        result_r_re[1][m0] += common_factor.x * ylm0_val * r_op_y;
                        result_r_im[1][m0] += common_factor.y * ylm0_val * r_op_y;
                        result_r_re[2][m0] += common_factor.x * ylm0_val * r_op_z;
                        result_r_im[2][m0] += common_factor.y * ylm0_val * r_op_z;
                    }
                }
            }
        } // End radial loop
    }     // End angular loop

    //-------------------------------------------------------------------------
    // Parallel reduction and output
    // Uses warp shuffle for efficiency, followed by cross-warp reduction
    //-------------------------------------------------------------------------

    const int out_base = norb_idx * nlm_dim * natomwfc;
    const int warp_id = tid / 32;
    const int lane_id = tid % 32;

    for (int m0 = 0; m0 < num_m0; m0++)
    {
        // Step 1: Warp-level reduction using shuffle
        double sum_re = warp_reduce_sum(result_re[m0]);
        double sum_im = warp_reduce_sum(result_im[m0]);

        // Step 2: First lane of each warp writes to shared memory
        if (lane_id == 0)
        {
            s_temp_re[warp_id] = sum_re;
            s_temp_im[warp_id] = sum_im;
        }
        __syncthreads();

        // Step 3: First warp reduces across all warps and writes output
        if (warp_id == 0)
        {
            sum_re = (lane_id < NUM_WARPS) ? s_temp_re[lane_id] : 0.0;
            sum_im = (lane_id < NUM_WARPS) ? s_temp_im[lane_id] : 0.0;
            sum_re = warp_reduce_sum(sum_re);
            sum_im = warp_reduce_sum(sum_im);

            if (lane_id == 0)
            {
                cuDoubleComplex result = make_cuDoubleComplex(sum_re, sum_im);
                result = cu_mul(result, exp_iAR0);
                result = cu_conj(result);
                nlm_out[out_base + 0 * natomwfc + m0_offset + m0] = result;
            }
        }
        __syncthreads();

        // Process current operator components (if nlm_dim == 4)
        if (nlm_dim == 4)
        {
            for (int d = 0; d < 3; d++)
            {
                double sum_r_re = warp_reduce_sum(result_r_re[d][m0]);
                double sum_r_im = warp_reduce_sum(result_r_im[d][m0]);

                if (lane_id == 0)
                {
                    s_temp_re[warp_id] = sum_r_re;
                    s_temp_im[warp_id] = sum_r_im;
                }
                __syncthreads();

                if (warp_id == 0)
                {
                    sum_r_re = (lane_id < NUM_WARPS) ? s_temp_re[lane_id] : 0.0;
                    sum_r_im = (lane_id < NUM_WARPS) ? s_temp_im[lane_id] : 0.0;
                    sum_r_re = warp_reduce_sum(sum_r_re);
                    sum_r_im = warp_reduce_sum(sum_r_im);

                    if (lane_id == 0)
                    {
                        cuDoubleComplex result_r = make_cuDoubleComplex(sum_r_re, sum_r_im);
                        result_r = cu_mul(result_r, exp_iAR0);
                        result_r = cu_conj(result_r);
                        nlm_out[out_base + (d + 1) * natomwfc + m0_offset + m0] = result_r;
                    }
                }
                __syncthreads();
            }
        }
    }
}

//=============================================================================
// Host-side Helper Functions
//=============================================================================

/**
 * @brief Copy integration grids to GPU constant memory
 *
 * Initializes the constant memory arrays with Lebedev-Laikov angular grid
 * and Gauss-Legendre radial grid for use in kernel integration.
 */
void copy_grids_to_device()
{
    // Copy Lebedev-Laikov 110-point angular quadrature grid
    CHECK_CUDA(cudaMemcpyToSymbol(d_lebedev_x,
                                  ModuleBase::Integral::Lebedev_Laikov_grid110_x,
                                  ANGULAR_GRID_NUM * sizeof(double)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_lebedev_y,
                                  ModuleBase::Integral::Lebedev_Laikov_grid110_y,
                                  ANGULAR_GRID_NUM * sizeof(double)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_lebedev_z,
                                  ModuleBase::Integral::Lebedev_Laikov_grid110_z,
                                  ANGULAR_GRID_NUM * sizeof(double)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_lebedev_w,
                                  ModuleBase::Integral::Lebedev_Laikov_grid110_w,
                                  ANGULAR_GRID_NUM * sizeof(double)));

    // Compute and copy Gauss-Legendre radial quadrature grid
    std::vector<double> h_gl_x(RADIAL_GRID_NUM);
    std::vector<double> h_gl_w(RADIAL_GRID_NUM);
    ModuleBase::Integral::Gauss_Legendre_grid_and_weight(RADIAL_GRID_NUM, h_gl_x.data(), h_gl_w.data());

    CHECK_CUDA(cudaMemcpyToSymbol(d_gl_x, h_gl_x.data(), RADIAL_GRID_NUM * sizeof(double)));
    CHECK_CUDA(cudaMemcpyToSymbol(d_gl_w, h_gl_w.data(), RADIAL_GRID_NUM * sizeof(double)));
}

} // namespace gpu
} // namespace module_rt
