#include "snap_psibeta_kernel.cuh"
#include "source_base/math_integral.h"

#include <cstdio>
#include <vector>

namespace module_rt
{
namespace gpu
{

//=============================================================================
// Constant Memory for Integration Grids
//=============================================================================

// Lebedev-Laikov angular grid (110 points)
__constant__ double d_lebedev_x[ANGULAR_GRID_NUM];
__constant__ double d_lebedev_y[ANGULAR_GRID_NUM];
__constant__ double d_lebedev_z[ANGULAR_GRID_NUM];
__constant__ double d_lebedev_w[ANGULAR_GRID_NUM];

// Gauss-Legendre radial grid (140 points)
__constant__ double d_gl_x[RADIAL_GRID_NUM];
__constant__ double d_gl_w[RADIAL_GRID_NUM];

//=============================================================================
// Spherical Harmonics Implementation
//=============================================================================

/**
 * @brief Compute real spherical harmonics Y_lm
 *        Based on the recursive formula used in ModuleBase::Ylm
 *        OPTIMIZED: Use atan2, sincos for better performance
 */
__device__ void compute_ylm_gpu(int L, double x, double y, double z, double* ylm)
{
    const double PI = 3.14159265358979323846;
    const double FOUR_PI = 4.0 * PI;
    const double SQRT2 = 1.41421356237309504880;

    // Y_00
    ylm[0] = 0.5 * sqrt(1.0 / PI);

    if (L == 0)
    {
        return;
    }

    // Compute r and cost
    double r2 = x * x + y * y + z * z;
    double r = sqrt(r2);

    double cost, sint, phi;
    if (r < 1e-10)
    {
        cost = 1.0;
        sint = 0.0;
        phi = 0.0;
    }
    else
    {
        cost = z / r;
        sint = sqrt(1.0 - cost * cost);
        // Use atan2 for robust phi computation (replaces multiple conditionals)
        phi = atan2(y, x);
    }

    // Ensure sint is non-negative (numerical safety)
    if (sint < 0.0)
    {
        sint = 0.0;
    }

    // Associated Legendre polynomials P_l^m
    double p[MAX_L + 1][MAX_L + 1] = {0};

    p[0][0] = 1.0;

    if (L >= 1)
    {
        p[1][0] = cost;
        p[1][1] = -sint;
    }

    for (int l = 2; l <= L; l++)
    {
        for (int m = 0; m <= l - 2; m++)
        {
            p[l][m] = ((2 * l - 1) * cost * p[l - 1][m] - (l - 1 + m) * p[l - 2][m]) / (double)(l - m);
        }
        p[l][l - 1] = (2 * l - 1) * cost * p[l - 1][l - 1];
        double fact = 1.0;
        for (int i = 1; i <= 2 * l - 1; i += 2)
        {
            fact *= i;
        }
        double sint_pow = 1.0;
        for (int i = 0; i < l; i++)
        {
            sint_pow *= sint;
        }
        p[l][l] = fact * sint_pow;
        if (l % 2 == 1)
        {
            p[l][l] = -p[l][l];
        }
    }

    // Compute Y_lm
    int lm = 0;
    for (int l = 0; l <= L; l++)
    {
        double c = sqrt((2.0 * l + 1.0) / FOUR_PI);

        ylm[lm] = c * p[l][0];
        lm++;

        for (int m = 1; m <= l; m++)
        {
            double fact_ratio = 1.0;
            for (int i = l - m + 1; i <= l + m; i++)
            {
                fact_ratio *= i;
            }
            double same = c * sqrt(1.0 / fact_ratio) * SQRT2;

            // Use sincos for efficiency (computes both sin and cos together)
            double sin_mphi, cos_mphi;
            sincos(m * phi, &sin_mphi, &cos_mphi);

            ylm[lm] = same * p[l][m] * cos_mphi;
            lm++;

            ylm[lm] = same * p[l][m] * sin_mphi;
            lm++;
        }
    }
}

//=============================================================================
// Warp Reduction for double
//=============================================================================

__device__ __forceinline__ double warp_reduce_sum(double val)
{
    for (int offset = 16; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(0xffffffff, val, offset);
    }
    return val;
}

//=============================================================================
// Atom-Level Batch Kernel Implementation
// Processes ALL neighbors for a center atom in a single kernel launch
//=============================================================================

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
    int norb_idx = blockIdx.x; // Which (neighbor, orbital) pair
    int proj_idx = blockIdx.y; // Which projector
    int tid = threadIdx.x;

    if (norb_idx >= total_neighbor_orbitals || proj_idx >= nproj)
    {
        return;
    }

    // Load neighbor-orbital data (includes R1)
    const NeighborOrbitalData& norb = neighbor_orbitals[norb_idx];
    const ProjectorData& proj = projectors[proj_idx];

    double3 R1 = norb.R1;
    int L1 = norb.L1;
    int m1 = norb.m1;
    int L0 = proj.L0;
    int m0_offset = proj_m0_offset[proj_idx];

    // Skip if L values exceed our precomputed limit
    if (L1 > MAX_L || L0 > MAX_L)
    {
        return;
    }

    // Compute geometry
    double3 dR = make_double3(R1.x - R0.x, R1.y - R0.y, R1.z - R0.z);
    double distance01 = sqrt(dR.x * dR.x + dR.y * dR.y + dR.z * dR.z);

    double r1_max = norb.psi_rcut;
    // Integration range: use projector's radial grid range [r_min, r_max]
    // This matches the CPU implementation in snap_psibeta_half_tddft.cpp
    double r_min = proj.r_min;
    double r_max = proj.r_max;
    double xl = 0.5 * (r_max - r_min);
    double xmean = 0.5 * (r_max + r_min);

    // Phase factor exp(i A · R0)
    double AdotR0 = A.x * R0.x + A.y * R0.y + A.z * R0.z;
    cuDoubleComplex exp_iAR0 = cu_exp_i(AdotR0);

    // Shared memory for reduction
    __shared__ double s_temp_re[BLOCK_SIZE];
    __shared__ double s_temp_im[BLOCK_SIZE];

    // Initialize accumulators in registers
    double result_re[MAX_YLM_SIZE];
    double result_im[MAX_YLM_SIZE];
    double result_r_re[3][MAX_YLM_SIZE];
    double result_r_im[3][MAX_YLM_SIZE];

    int num_m0 = 2 * L0 + 1;
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

    // Main integration loop - RESTRUCTURED: angular loop outside, radial inside
    // ylm0 only depends on angular direction (leb_x,y,z), not radial distance
    // This saves 140x redundant ylm0 computations per angular point
    for (int ian = tid; ian < ANGULAR_GRID_NUM; ian += BLOCK_SIZE)
    {
        double leb_x = d_lebedev_x[ian];
        double leb_y = d_lebedev_y[ian];
        double leb_z = d_lebedev_z[ian];
        double w_ang = d_lebedev_w[ian];

        // Precompute ylm0 ONCE per angular point (doesn't depend on r_val)
        double ylm0[MAX_YLM_SIZE];
        compute_ylm_gpu(L0, leb_x, leb_y, leb_z, ylm0);
        int offset_L0 = L0 * L0;

        // Precompute A dot direction (doesn't depend on r_val)
        double A_dot_leb = A.x * leb_x + A.y * leb_y + A.z * leb_z;

        // Precompute dR components
        double dRx = R0.x - R1.x;
        double dRy = R0.y - R1.y;
        double dRz = R0.z - R1.z;

#pragma unroll 4
        for (int ir = 0; ir < RADIAL_GRID_NUM; ir++)
        {
            double r_val = xmean + xl * d_gl_x[ir];
            double w_rad = xl * d_gl_w[ir];

            // Local position relative to R0
            double rx = r_val * leb_x;
            double ry = r_val * leb_y;
            double rz = r_val * leb_z;

            // Vector from R1 to integration point
            double tx = rx + dRx;
            double ty = ry + dRy;
            double tz = rz + dRz;
            double tnorm = sqrt(tx * tx + ty * ty + tz * tz);

            // Check psi cutoff
            if (tnorm <= r1_max)
            {
                // Compute Y_L1 (depends on direction to R1, varies with r_val)
                double ylm1[MAX_YLM_SIZE];
                if (tnorm > 1e-10)
                {
                    double inv_tnorm = 1.0 / tnorm;
                    compute_ylm_gpu(L1, tx * inv_tnorm, ty * inv_tnorm, tz * inv_tnorm, ylm1);
                }
                else
                {
                    compute_ylm_gpu(L1, 0.0, 0.0, 1.0, ylm1);
                }

                // Interpolate psi
                double psi_val
                    = interpolate_radial_gpu(psi_radial + norb.psi_offset, norb.psi_mesh, 1.0 / norb.psi_dk, tnorm);

                // Phase factor
                double phase = r_val * A_dot_leb;
                cuDoubleComplex exp_iAr = cu_exp_i(phase);

                // Interpolate beta
                double beta_val
                    = interpolate_radial_gpu(beta_radial + proj.beta_offset, proj.beta_mesh, 1.0 / proj.beta_dk, r_val);

                // Y_L1m1
                double ylm_L1_val = ylm1[L1 * L1 + m1];

                // Combined factor
                double factor = ylm_L1_val * psi_val * beta_val * r_val * w_rad * w_ang;
                cuDoubleComplex common_factor = cu_mul_real(exp_iAr, factor);

// Accumulate for all m0
#pragma unroll
                for (int m0 = 0; m0 < num_m0; m0++)
                {
                    double ylm0_val = ylm0[offset_L0 + m0];

                    result_re[m0] += common_factor.x * ylm0_val;
                    result_im[m0] += common_factor.y * ylm0_val;

                    if (nlm_dim == 4)
                    {
                        double r_op_x = rx + R0.x;
                        double r_op_y = ry + R0.y;
                        double r_op_z = rz + R0.z;

                        result_r_re[0][m0] += common_factor.x * ylm0_val * r_op_x;
                        result_r_im[0][m0] += common_factor.y * ylm0_val * r_op_x;
                        result_r_re[1][m0] += common_factor.x * ylm0_val * r_op_y;
                        result_r_im[1][m0] += common_factor.y * ylm0_val * r_op_y;
                        result_r_re[2][m0] += common_factor.x * ylm0_val * r_op_z;
                        result_r_im[2][m0] += common_factor.y * ylm0_val * r_op_z;
                    }
                }
            }
        }
    }

    // Reduction and output - OUTSIDE radial loop
    int out_base = norb_idx * nlm_dim * natomwfc;

    for (int m0 = 0; m0 < num_m0; m0++)
    {
        s_temp_re[tid] = result_re[m0];
        s_temp_im[tid] = result_im[m0];
        __syncthreads();

        for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1)
        {
            if (tid < s)
            {
                s_temp_re[tid] += s_temp_re[tid + s];
                s_temp_im[tid] += s_temp_im[tid + s];
            }
            __syncthreads();
        }

        if (tid == 0)
        {
            cuDoubleComplex result = make_cuDoubleComplex(s_temp_re[0], s_temp_im[0]);
            result = cu_mul(result, exp_iAR0);
            result = cu_conj(result);
            nlm_out[out_base + 0 * natomwfc + m0_offset + m0] = result;
        }
        __syncthreads();

        if (nlm_dim == 4)
        {
            for (int d = 0; d < 3; d++)
            {
                s_temp_re[tid] = result_r_re[d][m0];
                s_temp_im[tid] = result_r_im[d][m0];
                __syncthreads();

                for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1)
                {
                    if (tid < s)
                    {
                        s_temp_re[tid] += s_temp_re[tid + s];
                        s_temp_im[tid] += s_temp_im[tid + s];
                    }
                    __syncthreads();
                }

                if (tid == 0)
                {
                    cuDoubleComplex result_r = make_cuDoubleComplex(s_temp_re[0], s_temp_im[0]);
                    result_r = cu_mul(result_r, exp_iAR0);
                    result_r = cu_conj(result_r);
                    nlm_out[out_base + (d + 1) * natomwfc + m0_offset + m0] = result_r;
                }
                __syncthreads();
            }
        }
    }
}

//=============================================================================
// Host-side Helper Function Implementations
//=============================================================================

#define CUDA_CHECK_KERNEL(call)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        cudaError_t err = (call);                                                                                      \
        if (err != cudaSuccess)                                                                                        \
        {                                                                                                              \
            fprintf(stderr, "[CUDA Kernel] Error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));       \
        }                                                                                                              \
    } while (0)

void copy_grids_to_device()
{
    // Copy Lebedev-Laikov angular grid to constant memory
    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_lebedev_x,
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_x,
                                         ANGULAR_GRID_NUM * sizeof(double)));
    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_lebedev_y,
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_y,
                                         ANGULAR_GRID_NUM * sizeof(double)));
    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_lebedev_z,
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_z,
                                         ANGULAR_GRID_NUM * sizeof(double)));
    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_lebedev_w,
                                         ModuleBase::Integral::Lebedev_Laikov_grid110_w,
                                         ANGULAR_GRID_NUM * sizeof(double)));

    // Prepare and copy Gauss-Legendre radial grid to constant memory
    std::vector<double> h_gl_x(RADIAL_GRID_NUM);
    std::vector<double> h_gl_w(RADIAL_GRID_NUM);
    ModuleBase::Integral::Gauss_Legendre_grid_and_weight(RADIAL_GRID_NUM, h_gl_x.data(), h_gl_w.data());

    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_gl_x, h_gl_x.data(), RADIAL_GRID_NUM * sizeof(double)));
    CUDA_CHECK_KERNEL(cudaMemcpyToSymbol(d_gl_w, h_gl_w.data(), RADIAL_GRID_NUM * sizeof(double)));
}

} // namespace gpu
} // namespace module_rt
