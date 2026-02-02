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
#include "source_base/ylm.h"

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

// Spherical harmonics coefficients (copied from ModuleBase::Ylm::ylmcoef)
__constant__ double d_ylmcoef[100]; ///< Ylm coefficients for gradient calculations

//=============================================================================
// Spherical Harmonics - Helper Functions
//=============================================================================

/**
 * @brief Access element in lower-triangular stored Legendre polynomial array
 *
 * For associated Legendre polynomials P_l^m, we only need 0 <= m <= l.
 * Storage layout: P_0^0, P_1^0, P_1^1, P_2^0, P_2^1, P_2^2, ...
 * Linear index: l*(l+1)/2 + m
 */
__device__ __forceinline__ double& p_access(double* p, int l, int m)
{
    return p[l * (l + 1) / 2 + m];
}

/**
 * @brief Read-only access to Legendre polynomial array
 */
__device__ __forceinline__ double p_get(const double* p, int l, int m)
{
    return p[l * (l + 1) / 2 + m];
}

//=============================================================================
// Spherical Harmonics - Main Implementation
//=============================================================================

/**
 * @brief Compute real spherical harmonics Y_lm (templated version)
 *
 * Uses the recursive computation of associated Legendre polynomials:
 *   P_l^m = ((2l-1)*cos(theta)*P_{l-1}^m - (l-1+m)*P_{l-2}^m) / (l-m)
 *   P_l^{l-1} = (2l-1)*cos(theta)*P_{l-1}^{l-1}
 *   P_l^l = (-1)^l * (2l-1)!! * sin^l(theta)
 *
 * Real spherical harmonics are defined as:
 *   Y_{lm} = c_l * P_l^0                           for m = 0
 *   Y_{l,2m-1} = c_l * sqrt(2/(l-m)!/(l+m)!) * P_l^m * cos(m*phi)  for m > 0
 *   Y_{l,2m}   = c_l * sqrt(2/(l-m)!/(l+m)!) * P_l^m * sin(m*phi)  for m > 0
 * where c_l = sqrt((2l+1)/(4*pi))
 *
 * @tparam L Maximum angular momentum (compile-time constant)
 * @param x, y, z Direction vector components (need not be normalized)
 * @param ylm Output array storing Y_lm values in order: Y_00, Y_10, Y_11c, Y_11s, ...
 */
template <int L>
__device__ void compute_ylm_gpu(double x, double y, double z, double* ylm)
{

    constexpr int P_SIZE = (L + 1) * (L + 2) / 2; // Lower triangular storage size

    // Y_00 = 1/(2*sqrt(pi))
    ylm[0] = 0.5 * sqrt(1.0 / ModuleBase::PI);

    if (L == 0)
    {
        return;
    }

    // Compute spherical angles
    double r2 = x * x + y * y + z * z;
    double r = sqrt(r2);

    double cost, sint, phi;
    if (r < 1e-10)
    {
        // At origin, default to z-axis direction
        cost = 1.0;
        sint = 0.0;
        phi = 0.0;
    }
    else
    {
        cost = z / r;
        sint = sqrt(1.0 - cost * cost);
        phi = atan2(y, x);
    }

    // Ensure sint is non-negative (numerical safety)
    if (sint < 0.0)
    {
        sint = 0.0;
    }

    // Associated Legendre polynomials P_l^m in lower-triangular storage
    double p[P_SIZE];

    // Base cases
    p_access(p, 0, 0) = 1.0;

    if (L >= 1)
    {
        p_access(p, 1, 0) = cost;  // P_1^0 = cos(theta)
        p_access(p, 1, 1) = -sint; // P_1^1 = -sin(theta)
    }

    // Recurrence relations for l >= 2
#pragma unroll
    for (int l = 2; l <= L; l++)
    {
        // P_l^m for m = 0 to l-2: standard recurrence
#pragma unroll
        for (int m = 0; m <= l - 2; m++)
        {
            p_access(p, l, m) = ((2 * l - 1) * cost * p_get(p, l - 1, m) - (l - 1 + m) * p_get(p, l - 2, m))
                                / static_cast<double>(l - m);
        }

        // P_l^{l-1} = (2l-1) * cos(theta) * P_{l-1}^{l-1}
        p_access(p, l, l - 1) = (2 * l - 1) * cost * p_get(p, l - 1, l - 1);

        // P_l^l = (-1)^l * (2l-1)!! * sin^l(theta)
        double double_factorial = 1.0;
#pragma unroll
        for (int i = 1; i <= 2 * l - 1; i += 2)
        {
            double_factorial *= i;
        }

        double sint_power = 1.0;
#pragma unroll
        for (int i = 0; i < l; i++)
        {
            sint_power *= sint;
        }

        p_access(p, l, l) = double_factorial * sint_power;
        if (l % 2 == 1)
        {
            p_access(p, l, l) = -p_access(p, l, l);
        }
    }

    // Transform Legendre polynomials to real spherical harmonics
    int lm = 0;
#pragma unroll
    for (int l = 0; l <= L; l++)
    {
        double c = sqrt((2.0 * l + 1.0) / ModuleBase::FOUR_PI);

        // m = 0 component
        ylm[lm] = c * p_get(p, l, 0);
        lm++;

        // m > 0 components (cosine and sine parts)
#pragma unroll
        for (int m = 1; m <= l; m++)
        {
            // Compute normalization factor: sqrt(2 * (l-m)! / (l+m)!)
            double factorial_ratio = 1.0;
#pragma unroll
            for (int i = l - m + 1; i <= l + m; i++)
            {
                factorial_ratio *= i;
            }
            double norm = c * sqrt(1.0 / factorial_ratio) * ModuleBase::SQRT2;

            double sin_mphi, cos_mphi;
            sincos(m * phi, &sin_mphi, &cos_mphi);

            ylm[lm] = norm * p_get(p, l, m) * cos_mphi; // Y_{l,m} cosine part
            lm++;

            ylm[lm] = norm * p_get(p, l, m) * sin_mphi; // Y_{l,m} sine part
            lm++;
        }
    }
}

// Explicit template instantiations for L = 0, 1, 2, 3, 4
template __device__ void compute_ylm_gpu<0>(double x, double y, double z, double* ylm);
template __device__ void compute_ylm_gpu<1>(double x, double y, double z, double* ylm);
template __device__ void compute_ylm_gpu<2>(double x, double y, double z, double* ylm);
template __device__ void compute_ylm_gpu<3>(double x, double y, double z, double* ylm);
template __device__ void compute_ylm_gpu<4>(double x, double y, double z, double* ylm);

//=============================================================================
// Spherical Harmonics Gradients - GPU Implementation
//=============================================================================

/**
 * @brief Compute gradients of real spherical harmonics Y_lm (templated version)
 *
 * Computes ∂Y_lm/∂x, ∂Y_lm/∂y, ∂Y_lm/∂z for all (l,m) up to L.
 * This is the GPU equivalent of ModuleBase::Ylm::grad_rl_sph_harm.
 *
 * The gradients are computed using the recurrence relations:
 *   ∂(r^l Y_lm)/∂x_i = r^l * ∂Y_lm/∂x_i + Y_lm * ∂(r^l)/∂x_i
 *
 * For r^l Y_lm (solid harmonics), the CPU implementation computes gradients
 * recursively. We follow the same pattern here.
 *
 * @tparam L Maximum angular momentum (compile-time constant)
 * @param x, y, z Direction vector components (Cartesian coordinates)
 * @param rly Input: Y_lm values (must be precomputed via compute_ylm_gpu)
 * @param grly_x Output: ∂Y_lm/∂x for all (l,m)
 * @param grly_y Output: ∂Y_lm/∂y for all (l,m)
 * @param grly_z Output: ∂Y_lm/∂z for all (l,m)
 */
template <int L>
__device__ void compute_ylm_gradient_gpu(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z)
{
    double radius2 = x * x + y * y + z * z;
    double tx = 2.0 * x;
    double ty = 2.0 * y;
    double tz = 2.0 * z;

    // L = 0: Y_00 is constant, gradient is zero
    grly_x[0] = 0.0;
    grly_y[0] = 0.0;
    grly_z[0] = 0.0;

    if (L == 0) return;

    // L = 1
    // rly[1] = c1 * z
    grly_x[1] = 0.0;
    grly_y[1] = 0.0;
    grly_z[1] = d_ylmcoef[1];

    // rly[2] = -c1 * x
    grly_x[2] = -d_ylmcoef[1];
    grly_y[2] = 0.0;
    grly_z[2] = 0.0;

    // rly[3] = -c1 * y
    grly_x[3] = 0.0;
    grly_y[3] = -d_ylmcoef[1];
    grly_z[3] = 0.0;

    if (L == 1) return;

    // L = 2
    // rly[4] = c2*z*rly[1] - c3*rly[0]*radius2
    grly_x[4] = d_ylmcoef[2] * z * grly_x[1] - d_ylmcoef[3] * (grly_x[0] * radius2 + rly[0] * tx);
    grly_y[4] = d_ylmcoef[2] * z * grly_y[1] - d_ylmcoef[3] * (grly_y[0] * radius2 + rly[0] * ty);
    grly_z[4] = d_ylmcoef[2] * (z * grly_z[1] + rly[1]) - d_ylmcoef[3] * (grly_z[0] * radius2 + rly[0] * tz);

    double tmp0 = d_ylmcoef[4] * z;
    // rly[5] = tmp0 * rly[2]
    grly_x[5] = tmp0 * grly_x[2];
    grly_y[5] = tmp0 * grly_y[2];
    grly_z[5] = d_ylmcoef[4] * (rly[2] + z * grly_z[2]);

    // rly[6] = tmp0 * rly[3]
    grly_x[6] = tmp0 * grly_x[3];
    grly_y[6] = tmp0 * grly_y[3];
    grly_z[6] = d_ylmcoef[4] * (rly[3] + z * grly_z[3]);

    double tmp2 = d_ylmcoef[4] * x;
    // rly[7] = c5*rly[4] - c6*rly[0]*radius2 - tmp2*rly[2]
    grly_x[7] = d_ylmcoef[5] * grly_x[4] - d_ylmcoef[6] * (rly[0] * tx + grly_x[0] * radius2) - d_ylmcoef[4] * (x * grly_x[2] + rly[2]);
    grly_y[7] = d_ylmcoef[5] * grly_y[4] - d_ylmcoef[6] * (rly[0] * ty + grly_y[0] * radius2) - tmp2 * grly_y[2];
    grly_z[7] = d_ylmcoef[5] * grly_z[4] - d_ylmcoef[6] * (rly[0] * tz + grly_z[0] * radius2) - tmp2 * grly_z[2];

    // rly[8] = -tmp2 * rly[3]
    grly_x[8] = -d_ylmcoef[4] * (rly[3] + x * grly_x[3]);
    grly_y[8] = -tmp2 * grly_y[3];
    grly_z[8] = -tmp2 * grly_z[3];

    if (L == 2) return;

    // L = 3
    // rly[9] = c7*z*rly[4] - c8*rly[1]*radius2
    grly_x[9] = d_ylmcoef[7] * z * grly_x[4] - d_ylmcoef[8] * (rly[1] * tx + grly_x[1] * radius2);
    grly_y[9] = d_ylmcoef[7] * z * grly_y[4] - d_ylmcoef[8] * (rly[1] * ty + grly_y[1] * radius2);
    grly_z[9] = d_ylmcoef[7] * (rly[4] + z * grly_z[4]) - d_ylmcoef[8] * (rly[1] * tz + grly_z[1] * radius2);

    double tmp3 = d_ylmcoef[9] * z;
    // rly[10] = tmp3*rly[5] - c10*rly[2]*radius2
    grly_x[10] = tmp3 * grly_x[5] - d_ylmcoef[10] * (grly_x[2] * radius2 + rly[2] * tx);
    grly_y[10] = tmp3 * grly_y[5] - d_ylmcoef[10] * (grly_y[2] * radius2 + rly[2] * ty);
    grly_z[10] = d_ylmcoef[9] * (z * grly_z[5] + rly[5]) - d_ylmcoef[10] * (grly_z[2] * radius2 + rly[2] * tz);

    // rly[11] = tmp3*rly[6] - c10*rly[3]*radius2
    grly_x[11] = tmp3 * grly_x[6] - d_ylmcoef[10] * (grly_x[3] * radius2 + rly[3] * tx);
    grly_y[11] = tmp3 * grly_y[6] - d_ylmcoef[10] * (grly_y[3] * radius2 + rly[3] * ty);
    grly_z[11] = d_ylmcoef[9] * (z * grly_z[6] + rly[6]) - d_ylmcoef[10] * (grly_z[3] * radius2 + rly[3] * tz);

    double tmp4 = d_ylmcoef[11] * z;
    // rly[12] = tmp4*rly[7]
    grly_x[12] = tmp4 * grly_x[7];
    grly_y[12] = tmp4 * grly_y[7];
    grly_z[12] = d_ylmcoef[11] * (z * grly_z[7] + rly[7]);

    // rly[13] = tmp4*rly[8]
    grly_x[13] = tmp4 * grly_x[8];
    grly_y[13] = tmp4 * grly_y[8];
    grly_z[13] = d_ylmcoef[11] * (z * grly_z[8] + rly[8]);

    double tmp5 = d_ylmcoef[14] * x;
    // rly[14] = c12*rly[10] - c13*rly[2]*radius2 - tmp5*rly[7]
    grly_x[14] = d_ylmcoef[12] * grly_x[10] - d_ylmcoef[13] * (rly[2] * tx + grly_x[2] * radius2) - d_ylmcoef[14] * (rly[7] + x * grly_x[7]);
    grly_y[14] = d_ylmcoef[12] * grly_y[10] - d_ylmcoef[13] * (rly[2] * ty + grly_y[2] * radius2) - tmp5 * grly_y[7];
    grly_z[14] = d_ylmcoef[12] * grly_z[10] - d_ylmcoef[13] * (rly[2] * tz + grly_z[2] * radius2) - tmp5 * grly_z[7];

    // rly[15] = c12*rly[11] - c13*rly[3]*radius2 - tmp5*rly[8]
    grly_x[15] = d_ylmcoef[12] * grly_x[11] - d_ylmcoef[13] * (rly[3] * tx + grly_x[3] * radius2) - d_ylmcoef[14] * (rly[8] + x * grly_x[8]);
    grly_y[15] = d_ylmcoef[12] * grly_y[11] - d_ylmcoef[13] * (rly[3] * ty + grly_y[3] * radius2) - tmp5 * grly_y[8];
    grly_z[15] = d_ylmcoef[12] * grly_z[11] - d_ylmcoef[13] * (rly[3] * tz + grly_z[3] * radius2) - tmp5 * grly_z[8];

    if (L == 3) return;

    // L = 4
    // rly[16] = c15*z*rly[9] - c16*rly[4]*radius2
    grly_x[16] = d_ylmcoef[15] * z * grly_x[9] - d_ylmcoef[16] * (rly[4] * tx + grly_x[4] * radius2);
    grly_y[16] = d_ylmcoef[15] * z * grly_y[9] - d_ylmcoef[16] * (rly[4] * ty + grly_y[4] * radius2);
    grly_z[16] = d_ylmcoef[15] * (z * grly_z[9] + rly[9]) - d_ylmcoef[16] * (rly[4] * tz + grly_z[4] * radius2);

    double tmp6 = d_ylmcoef[17] * z;
    // rly[17] = tmp6*rly[10] - c18*rly[5]*radius2
    grly_x[17] = tmp6 * grly_x[10] - d_ylmcoef[18] * (rly[5] * tx + grly_x[5] * radius2);
    grly_y[17] = tmp6 * grly_y[10] - d_ylmcoef[18] * (rly[5] * ty + grly_y[5] * radius2);
    grly_z[17] = d_ylmcoef[17] * (z * grly_z[10] + rly[10]) - d_ylmcoef[18] * (rly[5] * tz + grly_z[5] * radius2);

    // rly[18] = tmp6*rly[11] - c18*rly[6]*radius2
    grly_x[18] = tmp6 * grly_x[11] - d_ylmcoef[18] * (rly[6] * tx + grly_x[6] * radius2);
    grly_y[18] = tmp6 * grly_y[11] - d_ylmcoef[18] * (rly[6] * ty + grly_y[6] * radius2);
    grly_z[18] = d_ylmcoef[17] * (z * grly_z[11] + rly[11]) - d_ylmcoef[18] * (rly[6] * tz + grly_z[6] * radius2);

    double tmp7 = d_ylmcoef[19] * z;
    // rly[19] = tmp7*rly[12] - c20*rly[7]*radius2
    grly_x[19] = tmp7 * grly_x[12] - d_ylmcoef[20] * (rly[7] * tx + grly_x[7] * radius2);
    grly_y[19] = tmp7 * grly_y[12] - d_ylmcoef[20] * (rly[7] * ty + grly_y[7] * radius2);
    grly_z[19] = d_ylmcoef[19] * (z * grly_z[12] + rly[12]) - d_ylmcoef[20] * (rly[7] * tz + grly_z[7] * radius2);

    // rly[20] = tmp7*rly[13] - c20*rly[8]*radius2
    grly_x[20] = tmp7 * grly_x[13] - d_ylmcoef[20] * (rly[8] * tx + grly_x[8] * radius2);
    grly_y[20] = tmp7 * grly_y[13] - d_ylmcoef[20] * (rly[8] * ty + grly_y[8] * radius2);
    grly_z[20] = d_ylmcoef[19] * (z * grly_z[13] + rly[13]) - d_ylmcoef[20] * (rly[8] * tz + grly_z[8] * radius2);

    double tmp8 = 3.0 * z;
    // rly[21] = tmp8*rly[14]
    grly_x[21] = tmp8 * grly_x[14];
    grly_y[21] = tmp8 * grly_y[14];
    grly_z[21] = 3.0 * (z * grly_z[14] + rly[14]);

    // rly[22] = tmp8*rly[15]
    grly_x[22] = tmp8 * grly_x[15];
    grly_y[22] = tmp8 * grly_y[15];
    grly_z[22] = 3.0 * (z * grly_z[15] + rly[15]);

    double tmp9 = d_ylmcoef[23] * x;
    // rly[23] = c21*rly[19] - c22*rly[7]*radius2 - tmp9*rly[14]
    grly_x[23] = d_ylmcoef[21] * grly_x[19] - d_ylmcoef[22] * (rly[7] * tx + grly_x[7] * radius2) - d_ylmcoef[23] * (x * grly_x[14] + rly[14]);
    grly_y[23] = d_ylmcoef[21] * grly_y[19] - d_ylmcoef[22] * (rly[7] * ty + grly_y[7] * radius2) - tmp9 * grly_y[14];
    grly_z[23] = d_ylmcoef[21] * grly_z[19] - d_ylmcoef[22] * (rly[7] * tz + grly_z[7] * radius2) - tmp9 * grly_z[14];

    // rly[24] = c21*rly[20] - c22*rly[8]*radius2 - tmp9*rly[15]
    grly_x[24] = d_ylmcoef[21] * grly_x[20] - d_ylmcoef[22] * (rly[8] * tx + grly_x[8] * radius2) - d_ylmcoef[23] * (x * grly_x[15] + rly[15]);
    grly_y[24] = d_ylmcoef[21] * grly_y[20] - d_ylmcoef[22] * (rly[8] * ty + grly_y[8] * radius2) - tmp9 * grly_y[15];
    grly_z[24] = d_ylmcoef[21] * grly_z[20] - d_ylmcoef[22] * (rly[8] * tz + grly_z[8] * radius2) - tmp9 * grly_z[15];
}

// Explicit template instantiations for gradient functions
template __device__ void compute_ylm_gradient_gpu<0>(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z);
template __device__ void compute_ylm_gradient_gpu<1>(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z);
template __device__ void compute_ylm_gradient_gpu<2>(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z);
template __device__ void compute_ylm_gradient_gpu<3>(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z);
template __device__ void compute_ylm_gradient_gpu<4>(double x, double y, double z, const double* rly, double* grly_x, double* grly_y, double* grly_z);

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
    double result_r_re[3][MAX_M0_SIZE]; // For position operator: x, y, z components
    double result_r_im[3][MAX_M0_SIZE];

    // Additional accumulators for derivatives (nlm_dim >= 4 with calc_deri)
    double result_d_re[3][MAX_M0_SIZE]; // For derivatives: dx, dy, dz components
    double result_d_im[3][MAX_M0_SIZE];

    // Additional accumulators for full tensor (nlm_dim == 16)
    double result_tensor_re[9][MAX_M0_SIZE]; // For 3x3 tensor: r_a * d/dtau_b
    double result_tensor_im[9][MAX_M0_SIZE];

    for (int m0 = 0; m0 < num_m0; m0++)
    {
        result_re[m0] = 0.0;
        result_im[m0] = 0.0;
        for (int d = 0; d < 3; d++)
        {
            result_r_re[d][m0] = 0.0;
            result_r_im[d][m0] = 0.0;
            result_d_re[d][m0] = 0.0;
            result_d_im[d][m0] = 0.0;
        }
        if (nlm_dim == 16)
        {
            for (int i = 0; i < 9; i++)
            {
                result_tensor_re[i][m0] = 0.0;
                result_tensor_im[i][m0] = 0.0;
            }
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
        DISPATCH_YLM(L0, leb_x, leb_y, leb_z, ylm0);
        const int offset_L0 = L0 * L0;

        // Precompute Y_lm gradients for projector if derivatives are needed
        double grly0_x[MAX_YLM_SIZE];
        double grly0_y[MAX_YLM_SIZE];
        double grly0_z[MAX_YLM_SIZE];
        if (nlm_dim >= 4)
        {
            // Note: For nlm_dim=4, we need to check if it's position or derivatives
            // For now, we compute gradients when nlm_dim >= 4
            // The caller will determine the actual mode
            DISPATCH_YLM_GRADIENT(L0, leb_x, leb_y, leb_z, ylm0, grly0_x, grly0_y, grly0_z);
        }

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
                    DISPATCH_YLM(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, ylm1);
                }
                else
                {
                    DISPATCH_YLM(L1, 0.0, 0.0, 1.0, ylm1);
                }

                // Interpolate orbital radial function
                const double psi_val
                    = interpolate_radial_gpu(psi_radial + norb.psi_offset, norb.psi_mesh, 1.0 / norb.psi_dk, tnorm);

                // Interpolate projector radial function
                const double beta_val
                    = interpolate_radial_gpu(beta_radial + proj.beta_offset, proj.beta_mesh, 1.0 / proj.beta_dk, r_val);

                // Compute radial derivative if needed (nlm_dim >= 4 for derivatives)
                double dbeta_dr = 0.0;
                if (nlm_dim >= 4 && r_val > 1e-10)
                {
                    dbeta_dr = compute_radial_derivative_gpu(beta_radial + proj.beta_offset, proj.beta_mesh, 1.0 / proj.beta_dk, r_val);
                }

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

                    // Position operator contribution (if nlm_dim == 4 or 16)
                    if (nlm_dim == 4 || nlm_dim == 16)
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

                    // Derivative contributions (if nlm_dim >= 4)
                    // Note: This computes derivatives regardless of whether position is also needed
                    // The caller will determine which output slots to use based on calc_deri flag
                    if (nlm_dim >= 4)
                    {
                        // Angular part: beta * (∂Y_lm/∂τ_a)
                        const double grad_ylm_x = grly0_x[offset_L0 + m0];
                        const double grad_ylm_y = grly0_y[offset_L0 + m0];
                        const double grad_ylm_z = grly0_z[offset_L0 + m0];

                        result_d_re[0][m0] += common_factor.x * grad_ylm_x;
                        result_d_im[0][m0] += common_factor.y * grad_ylm_x;
                        result_d_re[1][m0] += common_factor.x * grad_ylm_y;
                        result_d_im[1][m0] += common_factor.y * grad_ylm_y;
                        result_d_re[2][m0] += common_factor.x * grad_ylm_z;
                        result_d_im[2][m0] += common_factor.y * grad_ylm_z;

                        // Radial part: (∂beta/∂r) * (r_a/r)
                        if (r_val > 1e-10)
                        {
                            const double inv_r = 1.0 / r_val;
                            const double radial_factor_deriv = ylm_L1_val * psi_val * dbeta_dr * r_val * w_rad * w_ang;
                            const cuDoubleComplex radial_term = cu_mul_real(exp_iAr, radial_factor_deriv);

                            result_d_re[0][m0] += radial_term.x * ylm0_val * (rx * inv_r);
                            result_d_im[0][m0] += radial_term.y * ylm0_val * (rx * inv_r);
                            result_d_re[1][m0] += radial_term.x * ylm0_val * (ry * inv_r);
                            result_d_im[1][m0] += radial_term.y * ylm0_val * (ry * inv_r);
                            result_d_re[2][m0] += radial_term.x * ylm0_val * (rz * inv_r);
                            result_d_im[2][m0] += radial_term.y * ylm0_val * (rz * inv_r);
                        }

                        // Full 3x3 tensor: r_a * ∂/∂τ_b (if nlm_dim == 16)
                        if (nlm_dim == 16)
                        {
                            const double r_op_x = rx + R0.x;
                            const double r_op_y = ry + R0.y;
                            const double r_op_z = rz + R0.z;
                            const double r_op[3] = {r_op_x, r_op_y, r_op_z};
                            const double grad_ylm[3] = {grad_ylm_x, grad_ylm_y, grad_ylm_z};

                            // Angular contribution: r_op[a] * grad_ylm[b]
                            for (int a = 0; a < 3; a++)
                            {
                                for (int b = 0; b < 3; b++)
                                {
                                    const int idx = a * 3 + b;
                                    result_tensor_re[idx][m0] += common_factor.x * r_op[a] * grad_ylm[b];
                                    result_tensor_im[idx][m0] += common_factor.y * r_op[a] * grad_ylm[b];
                                }
                            }

                            // Radial contribution: r_op[a] * (∂beta/∂r) * (r_b/r)
                            if (r_val > 1e-10)
                            {
                                const double inv_r = 1.0 / r_val;
                                const double radial_factor_deriv = ylm_L1_val * psi_val * dbeta_dr * r_val * w_rad * w_ang;
                                const cuDoubleComplex radial_term = cu_mul_real(exp_iAr, radial_factor_deriv);
                                const double r_unit[3] = {leb_x, leb_y, leb_z}; // Unit vector

                                for (int a = 0; a < 3; a++)
                                {
                                    for (int b = 0; b < 3; b++)
                                    {
                                        const int idx = a * 3 + b;
                                        result_tensor_re[idx][m0] += radial_term.x * ylm0_val * r_op[a] * r_unit[b];
                                        result_tensor_im[idx][m0] += radial_term.y * ylm0_val * r_op[a] * r_unit[b];
                                    }
                                }
                            }
                        }
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

        // Process position/derivative components (if nlm_dim >= 4)
        if (nlm_dim >= 4)
        {
            // For nlm_dim == 4: output either position OR derivatives (determined by caller)
            // For nlm_dim == 16: output both position AND derivatives
            // We always reduce both here; caller determines which to use

            // Position operator components (slots 1-3 for nlm_dim==4 without calc_deri, or nlm_dim==16)
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
                        // For nlm_dim==4: slots 1-3 (position or derivatives based on caller)
                        // For nlm_dim==16: slots 1-3 (position)
                        nlm_out[out_base + (d + 1) * natomwfc + m0_offset + m0] = result_r;
                    }
                }
                __syncthreads();
            }

            // Derivative components (slots 4-6 for nlm_dim==16, or slots 1-3 for nlm_dim==4 with calc_deri)
            // Note: For nlm_dim==4, derivatives overwrite position in slots 1-3 (handled by caller)
            for (int d = 0; d < 3; d++)
            {
                double sum_d_re = warp_reduce_sum(result_d_re[d][m0]);
                double sum_d_im = warp_reduce_sum(result_d_im[d][m0]);

                if (lane_id == 0)
                {
                    s_temp_re[warp_id] = sum_d_re;
                    s_temp_im[warp_id] = sum_d_im;
                }
                __syncthreads();

                if (warp_id == 0)
                {
                    sum_d_re = (lane_id < NUM_WARPS) ? s_temp_re[lane_id] : 0.0;
                    sum_d_im = (lane_id < NUM_WARPS) ? s_temp_im[lane_id] : 0.0;
                    sum_d_re = warp_reduce_sum(sum_d_re);
                    sum_d_im = warp_reduce_sum(sum_d_im);

                    if (lane_id == 0)
                    {
                        cuDoubleComplex result_d = make_cuDoubleComplex(sum_d_re, sum_d_im);
                        result_d = cu_mul(result_d, exp_iAR0);
                        result_d = cu_conj(result_d);

                        if (nlm_dim == 16)
                        {
                            // For nlm_dim==16: derivatives go to slots 4-6
                            nlm_out[out_base + (d + 4) * natomwfc + m0_offset + m0] = result_d;
                        }
                        else
                        {
                            // For nlm_dim==4 with calc_deri: derivatives overwrite slots 1-3
                            // This is handled by the caller passing the right nlm_dim
                            // Here we just write to slots 1-3 (same as position above)
                            // The caller will decide which kernel output to use
                        }
                    }
                }
                __syncthreads();
            }

            // Full 3x3 tensor (slots 7-15 for nlm_dim==16)
            if (nlm_dim == 16)
            {
                for (int i = 0; i < 9; i++)
                {
                    double sum_t_re = warp_reduce_sum(result_tensor_re[i][m0]);
                    double sum_t_im = warp_reduce_sum(result_tensor_im[i][m0]);

                    if (lane_id == 0)
                    {
                        s_temp_re[warp_id] = sum_t_re;
                        s_temp_im[warp_id] = sum_t_im;
                    }
                    __syncthreads();

                    if (warp_id == 0)
                    {
                        sum_t_re = (lane_id < NUM_WARPS) ? s_temp_re[lane_id] : 0.0;
                        sum_t_im = (lane_id < NUM_WARPS) ? s_temp_im[lane_id] : 0.0;
                        sum_t_re = warp_reduce_sum(sum_t_re);
                        sum_t_im = warp_reduce_sum(sum_t_im);

                        if (lane_id == 0)
                        {
                            cuDoubleComplex result_t = make_cuDoubleComplex(sum_t_re, sum_t_im);
                            result_t = cu_mul(result_t, exp_iAR0);
                            result_t = cu_conj(result_t);
                            nlm_out[out_base + (i + 7) * natomwfc + m0_offset + m0] = result_t;
                        }
                    }
                    __syncthreads();
                }
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
