/*
 * NEP CUDA Compute - Core CUDA Kernels
 *
 * GPU implementation of NEP compute().
 * Three main kernels correspond to the three CPU compute steps:
 *
 *   1. nep_descriptor_kernel    - per-atom: descriptor + neural network → energy, Fp
 *   2. nep_force_radial_kernel  - per-neighbor-pair: radial force contribution
 *   3. nep_force_angular_kernel - per-neighbor-pair: angular force contribution
 *
 * These kernels operate on data already on the GPU.
 * They are designed to run after the neighbor list has been built on CPU
 * (neighbor list construction is complex and irregular on GPU; keeping it
 *  on CPU allows 2+3 to overlap with the next step's neighbor build via streams).
 */

#include "nep_cuda_compute.cuh"
#include <cuda_runtime.h>

// --------------- helper macros ---------------
#define CHECK_CUDA(call)                                      \
    do                                                        \
    {                                                         \
        cudaError_t e = (call);                               \
        if (e != cudaSuccess)                                 \
        {                                                     \
            fprintf(stderr, "CUDA error %s:%d: %s\n",        \
                    __FILE__, __LINE__, cudaGetErrorString(e));\
            exit(1);                                          \
        }                                                     \
    } while (0)

// =====================================================================
// Kernel 1: Descriptor + Neural Network (per-atom)
// =====================================================================

__global__ void nep_descriptor_kernel(
    // ---------- system ----------
    int N,                    // number of atoms
    int n_max_radial,         // paramb.n_max_radial
    int n_max_angular,        // paramb.n_max_angular
    int basis_size_radial,    // paramb.basis_size_radial
    int basis_size_angular,   // paramb.basis_size_angular
    int L_max,                // paramb.L_max
    int num_L,
    int num_types,            // paramb.num_types
    int num_types_sq,         // paramb.num_types_sq
    int num_c_radial,         // paramb.num_c_radial
    int dim,                  // annmb.dim
    int num_neurons1,         // annmb.num_neurons1
    int version,              // paramb.version
    // ---------- type ----------
    const int *g_type,
    // ---------- neighbor list ----------
    const int *g_NN_radial,
    const int *g_NL_radial,
    const int *g_NN_angular,
    const int *g_NL_angular,
    // ---------- pair distances ----------
    const double *g_x12_radial,
    const double *g_y12_radial,
    const double *g_z12_radial,
    const double *g_x12_angular,
    const double *g_y12_angular,
    const double *g_z12_angular,
    // ---------- radial cutoffs ----------
    const double *g_rc_radial,    // per-element [94]
    const double *g_rc_angular,   // per-element [94]
    // ---------- ANN parameters ----------
    const double *g_ann_c,        // expansion coefficients
    const double *g_w0,           // [annmb.num_para_ann] (all types packed)
    const double *g_b0,           // [annmb.num_para_ann]
    const double *g_w1,           // [annmb.num_para_ann]
    const double *g_b1,           // b1 pointer (size 1 or num_neurons1+1)
    // ---------- q scaler ----------
    const double *g_q_scaler,
    // ---------- output ----------
    double *g_potential,         // [N]
    double *g_Fp,                // [dim * N]  energy derivative wrt descriptor
    double *g_sum_fxyz            // [num_L * NUM_OF_ABC * N]
)
{
    int n1 = blockIdx.x * blockDim.x + threadIdx.x;
    if (n1 >= N) return;

    int t1 = g_type[n1];
    double q[NEP_CUDA_MAX_DIM] = {0.0};

    // ===== Part A: Radial descriptor =====
    int num_nn_radial = g_NN_radial[n1];
    for (int i1 = 0; i1 < num_nn_radial; ++i1)
    {
        int index = i1 * N + n1;
        int n2 = g_NL_radial[index];
        double r12[3] = {g_x12_radial[index], g_y12_radial[index], g_z12_radial[index]};
        double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

        int t2 = g_type[n2];
        double rc = (g_rc_radial[t1] + g_rc_radial[t2]) * 0.5;
        double rcinv = 1.0 / rc;

        double fc12;
        nep_cuda_find_fc(rc, rcinv, d12, fc12);
        double fn12[NEP_CUDA_MAX_NUM_N];
        nep_cuda_find_fn(basis_size_radial, rcinv, d12, fc12, fn12);

        for (int n = 0; n <= n_max_radial; ++n)
        {
            double gn12 = 0.0;
            for (int k = 0; k <= basis_size_radial; ++k)
            {
                int c_index = (n * (basis_size_radial + 1) + k) * num_types_sq;
                c_index += t1 * num_types + t2;
                gn12 += fn12[k] * g_ann_c[c_index];
            }
            q[n] += gn12;
        }
    }

    // ===== Part B: Angular descriptor =====
    int num_nn_angular = g_NN_angular[n1];
    for (int i1 = 0; i1 < num_nn_angular; ++i1)
    {
        int index = i1 * N + n1;
        int n2 = g_NL_angular[index];
        double r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
        double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

        int t2 = g_type[n2];
        double rc = (g_rc_angular[t1] + g_rc_angular[t2]) * 0.5;
        double rcinv = 1.0 / rc;

        double fc12;
        nep_cuda_find_fc(rc, rcinv, d12, fc12);
        double fn12[NEP_CUDA_MAX_NUM_N];
        nep_cuda_find_fn(basis_size_angular, rcinv, d12, fc12, fn12);

        // Per-order accumulation (n = 0..n_max_angular)
        for (int n = 0; n <= n_max_angular; ++n)
        {
            double gn12 = 0.0;
            for (int k = 0; k <= basis_size_angular; ++k)
            {
                int c_index = (n * (basis_size_angular + 1) + k) * num_types_sq;
                c_index += t1 * num_types + t2 + num_c_radial;
                gn12 += fn12[k] * g_ann_c[c_index];
            }
            // Accumulate S for this (n1, n2, order n)
            double s_order[NEP_CUDA_NUM_OF_ABC] = {0.0};
            nep_cuda_accumulate_s(L_max, d12, r12[0], r12[1], r12[2], gn12, s_order);

            // Q = f(S)
            nep_cuda_find_q(L_max, num_L, n_max_angular + 1, n, s_order, q + (n_max_radial + 1));

            // Save sum_fxyz for force kernels
            for (int abc = 0; abc < NEP_CUDA_NUM_OF_ABC; ++abc)
            {
                g_sum_fxyz[(n * NEP_CUDA_NUM_OF_ABC + abc) * N + n1] = s_order[abc];
            }
        }
    }

    // ===== Part C: Neural network (ANN) =====
    // Scale q by q_scaler
    for (int d = 0; d < dim; ++d)
    {
        q[d] = q[d] * g_q_scaler[d];
    }

    double F = 0.0;
    double Fp[NEP_CUDA_MAX_DIM] = {0.0};
    double latent_space[NEP_CUDA_MAX_NEURON] = {0.0};

    // Find weight pointers for this atom type
    // Weights are packed as: [num_types][num_neurons1 * dim] for w0/b0
    // w1 is [num_types][num_neurons1] or [num_neurons1] (shared across types)
    const double *w0_t1 = g_w0 + t1 * (num_neurons1 * dim);
    const double *b0_t1 = g_b0 + t1 * num_neurons1;

    if (version == 5)
    {
        nep_cuda_apply_ann_one_layer_nep5(
            dim, num_neurons1, w0_t1, b0_t1, g_w1, g_b1,
            q, F, Fp, latent_space);
    }
    else
    {
        nep_cuda_apply_ann_one_layer(
            dim, num_neurons1, w0_t1, b0_t1, g_w1,
            q, F, Fp, latent_space);
        F -= g_b1[0]; // subtract common bias for version < 5
    }

    g_potential[n1] = F;

    // Scale Fp by q_scaler for force computation
    for (int d = 0; d < dim; ++d)
    {
        g_Fp[d * N + n1] = Fp[d] * g_q_scaler[d];
    }
}

// =====================================================================
// Kernel 2: Radial Force (per-neighbor-pair)
// =====================================================================

__global__ void nep_force_radial_kernel(
    int N,
    int n_max_radial,
    int basis_size_radial,
    int num_types,
    int num_types_sq,
    int dim,
    int num_neurons1,
    int version,
    const int *g_type,
    const int *g_NN_radial,
    const int *g_NL_radial,
    const double *g_x12_radial,
    const double *g_y12_radial,
    const double *g_z12_radial,
    const double *g_rc_radial,
    const double *g_ann_c,
    const double *g_w0,
    const double *g_b0,
    const double *g_w1,
    const double *g_Fp,
    const double *g_q_scaler,
    double *g_fx,
    double *g_fy,
    double *g_fz,
    double *g_virial)
{
    // One thread per (atom, neighbor) pair
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N * NEP_CUDA_MN) return;

    int n1 = idx % N;
    int i1 = idx / N;
    if (i1 >= g_NN_radial[n1]) return;

    int t1 = g_type[n1];

    // Get neighbor info
    int index = i1 * N + n1;
    int n2 = g_NL_radial[index];
    int t2 = g_type[n2];
    double r12[3] = {g_x12_radial[index], g_y12_radial[index], g_z12_radial[index]};
    double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

    double rc = (g_rc_radial[t1] + g_rc_radial[t2]) * 0.5;
    double rcinv = 1.0 / rc;

    double fc12, fcp12;
    nep_cuda_find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);

    double fn12[NEP_CUDA_MAX_NUM_N], fnp12[NEP_CUDA_MAX_NUM_N];
    nep_cuda_find_fn_and_fnp(basis_size_radial, rcinv, d12, fc12, fcp12, fn12, fnp12);

    double d12inv = 1.0 / d12;

    // Load Fp for atom n1 (dim values)
    double Fp[NEP_CUDA_MAX_DIM];
    for (int d = 0; d < dim; ++d)
    {
        Fp[d] = g_Fp[d * N + n1];
    }

    // Accumulate force contribution from radial descriptor
    double fx = 0.0, fy = 0.0, fz = 0.0;
    double virial_xx = 0.0, virial_xy = 0.0, virial_xz = 0.0;
    double virial_yy = 0.0, virial_yz = 0.0, virial_zz = 0.0;

    for (int n = 0; n <= n_max_radial; ++n)
    {
        double Fp_n = Fp[n];
        if (Fp_n == 0.0) continue;

        double gnp12 = 0.0;
        for (int k = 0; k <= basis_size_radial; ++k)
        {
            int c_index = (n * (basis_size_radial + 1) + k) * num_types_sq;
            c_index += t1 * num_types + t2;
            gnp12 += fnp12[k] * g_ann_c[c_index];
        }

        double factor = Fp_n * gnp12 * d12inv;
        double dx = factor * r12[0];
        double dy = factor * r12[1];
        double dz = factor * r12[2];

        fx += dx;
        fy += dy;
        fz += dz;

        virial_xx += r12[0] * dx;
        virial_xy += r12[0] * dy;
        virial_xz += r12[0] * dz;
        virial_yy += r12[1] * dy;
        virial_yz += r12[1] * dz;
        virial_zz += r12[2] * dz;
    }

    // Write with atomicAdd (force contribution from pair)
    atomicAdd(&g_fx[n1], fx);
    atomicAdd(&g_fy[n1], fy);
    atomicAdd(&g_fz[n1], fz);
    atomicAdd(&g_fx[n2], -fx);
    atomicAdd(&g_fy[n2], -fy);
    atomicAdd(&g_fz[n2], -fz);

    // Virial (6 components)
    atomicAdd(&g_virial[n1], virial_xx);
    atomicAdd(&g_virial[n1 + N], virial_xy);
    atomicAdd(&g_virial[n1 + 2 * N], virial_xz);
    atomicAdd(&g_virial[n1 + 3 * N], virial_yy);
    atomicAdd(&g_virial[n1 + 4 * N], virial_yz);
    atomicAdd(&g_virial[n1 + 5 * N], virial_zz);

    // Symmetric virial contributions to n2
    atomicAdd(&g_virial[n2], virial_xx);
    atomicAdd(&g_virial[n2 + N], virial_xy);
    atomicAdd(&g_virial[n2 + 2 * N], virial_xz);
    atomicAdd(&g_virial[n2 + 3 * N], virial_yy);
    atomicAdd(&g_virial[n2 + 4 * N], virial_yz);
    atomicAdd(&g_virial[n2 + 5 * N], virial_zz);
}

// =====================================================================
// Kernel 3: Angular Force (per-neighbor-pair, same structure as radial)
// =====================================================================

__global__ void nep_force_angular_kernel(
    int N,
    int n_max_radial,
    int n_max_angular,
    int dim_angular,
    int basis_size_angular,
    int L_max,
    int num_L,
    int num_types,
    int num_types_sq,
    int num_c_radial,
    const int *g_type,
    const int *g_NN_angular,
    const int *g_NL_angular,
    const double *g_x12_angular,
    const double *g_y12_angular,
    const double *g_z12_angular,
    const double *g_rc_angular,
    const double *g_ann_c,
    const double *g_Fp,
    const double *g_sum_fxyz,
    double *g_fx,
    double *g_fy,
    double *g_fz,
    double *g_virial)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N * NEP_CUDA_MN) return;

    int n1 = idx % N;
    int i1 = idx / N;
    if (i1 >= g_NN_angular[n1]) return;

    int t1 = g_type[n1];
    int index = i1 * N + n1;
    int n2 = g_NL_angular[index];
    int t2 = g_type[n2];

    double r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
    double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);

    double rc = (g_rc_angular[t1] + g_rc_angular[t2]) * 0.5;
    double rcinv = 1.0 / rc;

    double fc12, fcp12;
    nep_cuda_find_fc_and_fcp(rc, rcinv, d12, fc12, fcp12);

    double fn12[NEP_CUDA_MAX_NUM_N], fnp12[NEP_CUDA_MAX_NUM_N];
    nep_cuda_find_fn_and_fnp(basis_size_angular, rcinv, d12, fc12, fcp12, fn12, fnp12);

    // Load Fp and sum_fxyz for atom n1 (per-CPU logic: load all orders)
    double Fp[NEP_CUDA_MAX_DIM_ANGULAR] = {0.0};
    double sum_fxyz[NEP_CUDA_NUM_OF_ABC * NEP_CUDA_MAX_NUM_N];
    for (int d = 0; d < dim_angular; ++d)
    {
        Fp[d] = g_Fp[(n_max_radial + 1 + d) * N + n1];
    }
    for (int d = 0; d < (n_max_angular + 1) * NEP_CUDA_NUM_OF_ABC; ++d)
    {
        sum_fxyz[d] = g_sum_fxyz[d * N + n1];
    }

    double f12[3] = {0.0};

    for (int n = 0; n <= n_max_angular; ++n)
    {
        double gn12 = 0.0;
        double gnp12 = 0.0;
        for (int k = 0; k <= basis_size_angular; ++k)
        {
            int c_index = (n * (basis_size_angular + 1) + k) * num_types_sq;
            c_index += t1 * num_types + t2 + num_c_radial;
            gn12 += fn12[k] * g_ann_c[c_index];
            gnp12 += fnp12[k] * g_ann_c[c_index];
        }

        // Full chain rule: d(q_angular)/d(x12) through S → Q → ANN
        nep_cuda_accumulate_f12(
            L_max, num_L, n, n_max_angular + 1, d12, r12, gn12, gnp12,
            Fp, sum_fxyz, f12);
    }

    // Write force contributions (Newton's third law)
    atomicAdd(&g_fx[n1], f12[0]);
    atomicAdd(&g_fy[n1], f12[1]);
    atomicAdd(&g_fz[n1], f12[2]);
    atomicAdd(&g_fx[n2], -f12[0]);
    atomicAdd(&g_fy[n2], -f12[1]);
    atomicAdd(&g_fz[n2], -f12[2]);

    // Virial contributions (matching CPU: g_virial[n2 + d * N] -= r12[d] * f12[d'])
    atomicAdd(&g_virial[n2 + 0 * N], -r12[0] * f12[0]);
    atomicAdd(&g_virial[n2 + 1 * N], -r12[0] * f12[1]);
    atomicAdd(&g_virial[n2 + 2 * N], -r12[0] * f12[2]);
    atomicAdd(&g_virial[n2 + 3 * N], -r12[1] * f12[0]);
    atomicAdd(&g_virial[n2 + 4 * N], -r12[1] * f12[1]);
    atomicAdd(&g_virial[n2 + 5 * N], -r12[1] * f12[2]);
    atomicAdd(&g_virial[n2 + 6 * N], -r12[2] * f12[0]);
    atomicAdd(&g_virial[n2 + 7 * N], -r12[2] * f12[1]);
    atomicAdd(&g_virial[n2 + 8 * N], -r12[2] * f12[2]);
}

// =====================================================================
// Kernel 4: ZBL Repulsive Force (per-neighbor-pair)
// =====================================================================

__global__ void nep_force_ZBL_kernel(
    int N,
    int num_types_zbl,
    int zbl_enabled,
    int zbl_flexible,
    double rc_inner,
    double rc_outer,
    int use_typewise_cutoff_zbl,
    double typewise_cutoff_zbl_factor,
    const int *g_type,
    const int *g_atomic_numbers,   // [num_types]
    const double *g_zbl_para,      // [num_types_sq_zbl * 10] or nullptr
    const int *g_NN_angular,
    const int *g_NL_angular,
    const double *g_x12_angular,
    const double *g_y12_angular,
    const double *g_z12_angular,
    double *g_fx,
    double *g_fy,
    double *g_fz,
    double *g_virial,
    double *g_pe) // potential energy per atom
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N * NEP_CUDA_MN) return;

    int n1 = idx % N;
    int i1 = idx / N;
    if (i1 >= g_NN_angular[n1]) return;

    int type1 = g_type[n1];
    int zi = g_atomic_numbers[type1] + 1;
    double pow_zi = pow(double(zi), 0.23);

    int index = i1 * N + n1;
    int n2 = g_NL_angular[index];
    double r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
    double d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
    double d12inv = 1.0 / d12;

    int type2 = g_type[n2];
    int zj = g_atomic_numbers[type2] + 1;
    double a_inv = (pow_zi + pow(double(zj), 0.23)) * 2.134563;
    double zizj = NEP_CUDA_K_C_SP * zi * zj;

    double f, fp;

    if (zbl_flexible)
    {
        int t1, t2;
        if (type1 < type2) { t1 = type1; t2 = type2; }
        else               { t1 = type2; t2 = type1; }
        int zbl_index = t1 * num_types_zbl - (t1 * (t1 - 1)) / 2 + (t2 - t1);
        double ZBL_para[10];
        for (int i = 0; i < 10; ++i)
        {
            ZBL_para[i] = g_zbl_para[10 * zbl_index + i];
        }
        nep_cuda_find_f_and_fp_zbl_flexible(ZBL_para, zizj, a_inv, d12, d12inv, f, fp);
    }
    else
    {
        double rc_i = rc_inner;
        double rc_o = rc_outer;
        if (use_typewise_cutoff_zbl)
        {
            rc_o = min(
                (nep_cuda_COVALENT_RADIUS[zi - 1] + nep_cuda_COVALENT_RADIUS[zj - 1]) * typewise_cutoff_zbl_factor,
                rc_o);
            rc_i = 0.0;
        }
        nep_cuda_find_f_and_fp_zbl(zizj, a_inv, rc_i, rc_o, d12, d12inv, f, fp);
    }

    double f2 = fp * d12inv * 0.5;
    double f12[3] = {r12[0] * f2, r12[1] * f2, r12[2] * f2};

    atomicAdd(&g_fx[n1], f12[0]);
    atomicAdd(&g_fy[n1], f12[1]);
    atomicAdd(&g_fz[n1], f12[2]);
    atomicAdd(&g_fx[n2], -f12[0]);
    atomicAdd(&g_fy[n2], -f12[1]);
    atomicAdd(&g_fz[n2], -f12[2]);

    atomicAdd(&g_virial[n2 + 0 * N], -r12[0] * f12[0]);
    atomicAdd(&g_virial[n2 + 1 * N], -r12[0] * f12[1]);
    atomicAdd(&g_virial[n2 + 2 * N], -r12[0] * f12[2]);
    atomicAdd(&g_virial[n2 + 3 * N], -r12[1] * f12[0]);
    atomicAdd(&g_virial[n2 + 4 * N], -r12[1] * f12[1]);
    atomicAdd(&g_virial[n2 + 5 * N], -r12[1] * f12[2]);
    atomicAdd(&g_virial[n2 + 6 * N], -r12[2] * f12[0]);
    atomicAdd(&g_virial[n2 + 7 * N], -r12[2] * f12[1]);
    atomicAdd(&g_virial[n2 + 8 * N], -r12[2] * f12[2]);

    if (g_pe)
    {
        atomicAdd(&g_pe[n1], f * 0.5);
    }
}

// =====================================================================
// Host-side workspace for persistent GPU buffers
// =====================================================================

struct NepCudaComputeWorkspace
{
    // Constant parameters (copy once, reuse across steps)
    double *d_rc_radial = nullptr;
    double *d_rc_angular = nullptr;
    double *d_ann_c = nullptr;
    double *d_w0 = nullptr;
    double *d_b0 = nullptr;
    double *d_w1 = nullptr;
    double *d_b1 = nullptr;
    double *d_q_scaler = nullptr;
    int *d_type = nullptr;

    // Variable data (per-step)
    int *d_NN_radial = nullptr;
    int *d_NL_radial = nullptr;
    int *d_NN_angular = nullptr;
    int *d_NL_angular = nullptr;
    double *d_x12_radial = nullptr;
    double *d_y12_radial = nullptr;
    double *d_z12_radial = nullptr;
    double *d_x12_angular = nullptr;
    double *d_y12_angular = nullptr;
    double *d_z12_angular = nullptr;

    // Output-intermediate
    double *d_potential = nullptr;
    double *d_Fp = nullptr;
    double *d_sum_fxyz = nullptr;

    // Force output
    double *d_fx = nullptr;
    double *d_fy = nullptr;
    double *d_fz = nullptr;
    double *d_virial = nullptr;

    int capacity = 0;
    bool params_loaded = false;
};

// =====================================================================
// Timing breakdown structure
// =====================================================================

struct NepCudaComputeTiming
{
    float h2d_copy_ms = 0.0f;      // Host→Device data transfer
    float descriptor_ms = 0.0f;    // Kernel 1: descriptor + ANN
    float force_radial_ms = 0.0f;  // Kernel 2: radial force
    float force_angular_ms = 0.0f; // Kernel 3: angular force
    float d2h_copy_ms = 0.0f;      // Device→Host result transfer
    float total_ms = 0.0f;         // Total GPU time
};

static void time_event_ms(cudaEvent_t start, cudaEvent_t stop, float &ms)
{
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&ms, start, stop);
}

// =====================================================================
// Main GPU compute entry point (untimed version)
// =====================================================================

void nep_cuda_compute(
    int N,
    const int *type,
    const int *NN_radial, const int *NL_radial,
    const int *NN_angular, const int *NL_angular,
    const double *x12_radial, const double *y12_radial, const double *z12_radial,
    const double *x12_angular, const double *y12_angular, const double *z12_angular,
    // NEP parameters
    int n_max_radial, int n_max_angular,
    int basis_size_radial, int basis_size_angular,
    int L_max, int num_L, int num_types, int num_types_sq, int num_c_radial,
    int dim, int num_neurons1, int version,
    const double *rc_radial, const double *rc_angular,
    const double *ann_c, int num_para,
    const double *w0, const double *b0, const double *w1, const double *b1,
    const double *q_scaler,
    // Output
    double *potential, double *force, double *virial)
{
    int size_type = N * sizeof(int);
    int size_N = N * sizeof(int);
    int size_double_N = N * sizeof(double);

    int MN = NEP_CUDA_MN;
    int size_nl = N * MN * sizeof(int);
    int size_nl_d = N * MN * sizeof(double);

    int *d_type, *d_NN_r, *d_NL_r, *d_NN_a, *d_NL_a;
    double *d_x12_r, *d_y12_r, *d_z12_r;
    double *d_x12_a, *d_y12_a, *d_z12_a;
    double *d_rc_r, *d_rc_a, *d_ann_c, *d_w0, *d_b0, *d_w1, *d_b1, *d_qs;
    double *d_pot, *d_Fp, *d_sfxyz, *d_fx, *d_fy, *d_fz, *d_vir;

    CHECK_CUDA(cudaMalloc(&d_type, size_type));
    CHECK_CUDA(cudaMalloc(&d_NN_r, size_N));
    CHECK_CUDA(cudaMalloc(&d_NL_r, size_nl));
    CHECK_CUDA(cudaMalloc(&d_NN_a, size_N));
    CHECK_CUDA(cudaMalloc(&d_NL_a, size_nl));
    CHECK_CUDA(cudaMalloc(&d_x12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_y12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_z12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_x12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_y12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_z12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_rc_r, 94 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rc_a, 94 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_ann_c, num_para * sizeof(double)));
    int w_size = num_types * num_neurons1 * dim * sizeof(double);
    int b_size = num_types * num_neurons1 * sizeof(double);
    int w1_size = num_types * num_neurons1 * sizeof(double);
    CHECK_CUDA(cudaMalloc(&d_w0, w_size));
    CHECK_CUDA(cudaMalloc(&d_b0, b_size));
    CHECK_CUDA(cudaMalloc(&d_w1, w1_size));
    CHECK_CUDA(cudaMalloc(&d_b1, (num_neurons1 + 1) * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_qs, dim * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_pot, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_Fp, dim * N * sizeof(double)));
    int sfxyz_size = num_L * NEP_CUDA_NUM_OF_ABC * N * sizeof(double);
    CHECK_CUDA(cudaMalloc(&d_sfxyz, sfxyz_size));
    CHECK_CUDA(cudaMalloc(&d_fx, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_fy, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_fz, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_vir, 9 * N * sizeof(double)));

    // Copy input data H2D
    CHECK_CUDA(cudaMemcpy(d_type, type, size_type, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NN_r, NN_radial, size_N, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NL_r, NL_radial, size_nl, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NN_a, NN_angular, size_N, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NL_a, NL_angular, size_nl, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_x12_r, x12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_y12_r, y12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_z12_r, z12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_x12_a, x12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_y12_a, y12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_z12_a, z12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rc_r, rc_radial, 94 * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rc_a, rc_angular, 94 * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_ann_c, ann_c, num_para * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_w0, w0, w_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_b0, b0, b_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_w1, w1, w1_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_b1, b1, (num_neurons1 + 1) * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_qs, q_scaler, dim * sizeof(double), cudaMemcpyHostToDevice));

    // Zero output buffers
    CHECK_CUDA(cudaMemset(d_pot, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_Fp, 0, dim * N * sizeof(double)));
    CHECK_CUDA(cudaMemset(d_sfxyz, 0, sfxyz_size));
    CHECK_CUDA(cudaMemset(d_fx, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_fy, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_fz, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_vir, 0, 9 * N * sizeof(double)));

    // ---- Launch kernel 1: descriptor ----
    int block_size = 128;
    int grid_size = (N + block_size - 1) / block_size;
    nep_descriptor_kernel<<<grid_size, block_size>>>(
        N, n_max_radial, n_max_angular,
        basis_size_radial, basis_size_angular,
        L_max, num_L, num_types, num_types_sq, num_c_radial,
        dim, num_neurons1, version,
        d_type,
        d_NN_r, d_NL_r, d_NN_a, d_NL_a,
        d_x12_r, d_y12_r, d_z12_r,
        d_x12_a, d_y12_a, d_z12_a,
        d_rc_r, d_rc_a,
        d_ann_c, d_w0, d_b0, d_w1, d_b1,
        d_qs,
        d_pot, d_Fp, d_sfxyz);
    cudaDeviceSynchronize();

    // ---- Launch kernel 2: radial force ----
    int total_pairs_radial = N * MN;
    grid_size = (total_pairs_radial + block_size - 1) / block_size;
    nep_force_radial_kernel<<<grid_size, block_size>>>(
        N, n_max_radial, basis_size_radial,
        num_types, num_types_sq, dim, num_neurons1, version,
        d_type, d_NN_r, d_NL_r,
        d_x12_r, d_y12_r, d_z12_r,
        d_rc_r, d_ann_c,
        d_w0, d_b0, d_w1,
        d_Fp, d_qs,
        d_fx, d_fy, d_fz, d_vir);
    cudaDeviceSynchronize();

    // ---- Launch kernel 3: angular force ----
    int total_pairs_angular = N * MN;
    grid_size = (total_pairs_angular + block_size - 1) / block_size;
    nep_force_angular_kernel<<<grid_size, block_size>>>(
        N, n_max_angular, basis_size_angular,
        L_max, num_L, num_types, num_types_sq, num_c_radial, dim,
        d_type, d_NN_a, d_NL_a,
        d_x12_a, d_y12_a, d_z12_a,
        d_rc_a, d_ann_c,
        d_Fp, d_sfxyz,
        d_fx, d_fy, d_fz, d_vir);
    cudaDeviceSynchronize();

    // Copy results D2H
    CHECK_CUDA(cudaMemcpy(potential, d_pot, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force, d_fx, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force + N, d_fy, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force + 2 * N, d_fz, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(virial, d_vir, 9 * N * sizeof(double), cudaMemcpyDeviceToHost));

    // Cleanup
    cudaFree(d_type);
    cudaFree(d_NN_r); cudaFree(d_NL_r);
    cudaFree(d_NN_a); cudaFree(d_NL_a);
    cudaFree(d_x12_r); cudaFree(d_y12_r); cudaFree(d_z12_r);
    cudaFree(d_x12_a); cudaFree(d_y12_a); cudaFree(d_z12_a);
    cudaFree(d_rc_r); cudaFree(d_rc_a);
    cudaFree(d_ann_c);
    cudaFree(d_w0); cudaFree(d_b0); cudaFree(d_w1); cudaFree(d_b1);
    cudaFree(d_qs);
    cudaFree(d_pot); cudaFree(d_Fp); cudaFree(d_sfxyz);
    cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz); cudaFree(d_vir);
}

// =====================================================================
// Timed version with CUDA Event profiling
// =====================================================================

void nep_cuda_compute_timed(
    int N,
    const int *type,
    const int *NN_radial, const int *NL_radial,
    const int *NN_angular, const int *NL_angular,
    const double *x12_radial, const double *y12_radial, const double *z12_radial,
    const double *x12_angular, const double *y12_angular, const double *z12_angular,
    int n_max_radial, int n_max_angular,
    int basis_size_radial, int basis_size_angular,
    int L_max, int num_L, int num_types, int num_types_sq, int num_c_radial,
    int dim, int num_neurons1, int version,
    const double *rc_radial, const double *rc_angular,
    const double *ann_c, int num_para,
    const double *w0, const double *b0, const double *w1, const double *b1,
    const double *q_scaler,
    double *potential, double *force, double *virial,
    NepCudaComputeTiming &timing)
{
    // Create CUDA events for timing
    cudaEvent_t ev_total_start, ev_total_stop;
    cudaEvent_t ev_h2d_start, ev_h2d_stop;
    cudaEvent_t ev_desc_start, ev_desc_stop;
    cudaEvent_t ev_fr_start, ev_fr_stop;
    cudaEvent_t ev_fa_start, ev_fa_stop;
    cudaEvent_t ev_d2h_start, ev_d2h_stop;

    cudaEventCreate(&ev_total_start);
    cudaEventCreate(&ev_total_stop);
    cudaEventCreate(&ev_h2d_start);
    cudaEventCreate(&ev_h2d_stop);
    cudaEventCreate(&ev_desc_start);
    cudaEventCreate(&ev_desc_stop);
    cudaEventCreate(&ev_fr_start);
    cudaEventCreate(&ev_fr_stop);
    cudaEventCreate(&ev_fa_start);
    cudaEventCreate(&ev_fa_stop);
    cudaEventCreate(&ev_d2h_start);
    cudaEventCreate(&ev_d2h_stop);

    cudaEventRecord(ev_total_start);

    int size_type = N * sizeof(int);
    int size_N = N * sizeof(int);
    int size_double_N = N * sizeof(double);
    int MN = NEP_CUDA_MN;
    int size_nl = N * MN * sizeof(int);
    int size_nl_d = N * MN * sizeof(double);

    int *d_type, *d_NN_r, *d_NL_r, *d_NN_a, *d_NL_a;
    double *d_x12_r, *d_y12_r, *d_z12_r;
    double *d_x12_a, *d_y12_a, *d_z12_a;
    double *d_rc_r, *d_rc_a, *d_ann_c, *d_w0, *d_b0, *d_w1, *d_b1, *d_qs;
    double *d_pot, *d_Fp, *d_sfxyz, *d_fx, *d_fy, *d_fz, *d_vir;

    CHECK_CUDA(cudaMalloc(&d_type, size_type));
    CHECK_CUDA(cudaMalloc(&d_NN_r, size_N));
    CHECK_CUDA(cudaMalloc(&d_NL_r, size_nl));
    CHECK_CUDA(cudaMalloc(&d_NN_a, size_N));
    CHECK_CUDA(cudaMalloc(&d_NL_a, size_nl));
    CHECK_CUDA(cudaMalloc(&d_x12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_y12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_z12_r, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_x12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_y12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_z12_a, size_nl_d));
    CHECK_CUDA(cudaMalloc(&d_rc_r, 94 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rc_a, 94 * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_ann_c, num_para * sizeof(double)));
    int w_size = num_types * num_neurons1 * dim * sizeof(double);
    int b_size = num_types * num_neurons1 * sizeof(double);
    int w1_size = num_types * num_neurons1 * sizeof(double);
    CHECK_CUDA(cudaMalloc(&d_w0, w_size));
    CHECK_CUDA(cudaMalloc(&d_b0, b_size));
    CHECK_CUDA(cudaMalloc(&d_w1, w1_size));
    CHECK_CUDA(cudaMalloc(&d_b1, (num_neurons1 + 1) * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_qs, dim * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_pot, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_Fp, dim * N * sizeof(double)));
    int sfxyz_size = num_L * NEP_CUDA_NUM_OF_ABC * N * sizeof(double);
    CHECK_CUDA(cudaMalloc(&d_sfxyz, sfxyz_size));
    CHECK_CUDA(cudaMalloc(&d_fx, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_fy, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_fz, size_double_N));
    CHECK_CUDA(cudaMalloc(&d_vir, 9 * N * sizeof(double)));

    // === Phase 1: H2D copy ===
    cudaEventRecord(ev_h2d_start);
    CHECK_CUDA(cudaMemcpy(d_type, type, size_type, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NN_r, NN_radial, size_N, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NL_r, NL_radial, size_nl, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NN_a, NN_angular, size_N, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_NL_a, NL_angular, size_nl, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_x12_r, x12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_y12_r, y12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_z12_r, z12_radial, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_x12_a, x12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_y12_a, y12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_z12_a, z12_angular, size_nl_d, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rc_r, rc_radial, 94 * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_rc_a, rc_angular, 94 * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_ann_c, ann_c, num_para * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_w0, w0, w_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_b0, b0, b_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_w1, w1, w1_size, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_b1, b1, (num_neurons1 + 1) * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_qs, q_scaler, dim * sizeof(double), cudaMemcpyHostToDevice));
    cudaEventRecord(ev_h2d_stop);

    // Zero output buffers
    CHECK_CUDA(cudaMemset(d_pot, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_Fp, 0, dim * N * sizeof(double)));
    CHECK_CUDA(cudaMemset(d_sfxyz, 0, sfxyz_size));
    CHECK_CUDA(cudaMemset(d_fx, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_fy, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_fz, 0, size_double_N));
    CHECK_CUDA(cudaMemset(d_vir, 0, 9 * N * sizeof(double)));

    int block_size = 128;

    // === Phase 2: Descriptor + ANN kernel ===
    cudaEventRecord(ev_desc_start);
    int grid_size = (N + block_size - 1) / block_size;
    nep_descriptor_kernel<<<grid_size, block_size>>>(
        N, n_max_radial, n_max_angular,
        basis_size_radial, basis_size_angular,
        L_max, num_L, num_types, num_types_sq, num_c_radial,
        dim, num_neurons1, version,
        d_type,
        d_NN_r, d_NL_r, d_NN_a, d_NL_a,
        d_x12_r, d_y12_r, d_z12_r,
        d_x12_a, d_y12_a, d_z12_a,
        d_rc_r, d_rc_a,
        d_ann_c, d_w0, d_b0, d_w1, d_b1,
        d_qs,
        d_pot, d_Fp, d_sfxyz);
    cudaDeviceSynchronize();
    cudaEventRecord(ev_desc_stop);

    // === Phase 3: Radial force kernel ===
    cudaEventRecord(ev_fr_start);
    int total_pairs_radial = N * MN;
    grid_size = (total_pairs_radial + block_size - 1) / block_size;
    nep_force_radial_kernel<<<grid_size, block_size>>>(
        N, n_max_radial, basis_size_radial,
        num_types, num_types_sq, dim, num_neurons1, version,
        d_type, d_NN_r, d_NL_r,
        d_x12_r, d_y12_r, d_z12_r,
        d_rc_r, d_ann_c,
        d_w0, d_b0, d_w1,
        d_Fp, d_qs,
        d_fx, d_fy, d_fz, d_vir);
    cudaDeviceSynchronize();
    cudaEventRecord(ev_fr_stop);

    // === Phase 4: Angular force kernel ===
    cudaEventRecord(ev_fa_start);
    int total_pairs_angular = N * MN;
    grid_size = (total_pairs_angular + block_size - 1) / block_size;
    nep_force_angular_kernel<<<grid_size, block_size>>>(
        N, n_max_angular, basis_size_angular,
        L_max, num_L, num_types, num_types_sq, num_c_radial, dim,
        d_type, d_NN_a, d_NL_a,
        d_x12_a, d_y12_a, d_z12_a,
        d_rc_a, d_ann_c,
        d_Fp, d_sfxyz,
        d_fx, d_fy, d_fz, d_vir);
    cudaDeviceSynchronize();
    cudaEventRecord(ev_fa_stop);

    // === Phase 5: D2H copy ===
    cudaEventRecord(ev_d2h_start);
    CHECK_CUDA(cudaMemcpy(potential, d_pot, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force, d_fx, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force + N, d_fy, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(force + 2 * N, d_fz, size_double_N, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(virial, d_vir, 9 * N * sizeof(double), cudaMemcpyDeviceToHost));
    cudaEventRecord(ev_d2h_stop);

    cudaEventRecord(ev_total_stop);

    // Extract timing
    time_event_ms(ev_h2d_start, ev_h2d_stop, timing.h2d_copy_ms);
    time_event_ms(ev_desc_start, ev_desc_stop, timing.descriptor_ms);
    time_event_ms(ev_fr_start, ev_fr_stop, timing.force_radial_ms);
    time_event_ms(ev_fa_start, ev_fa_stop, timing.force_angular_ms);
    time_event_ms(ev_d2h_start, ev_d2h_stop, timing.d2h_copy_ms);
    time_event_ms(ev_total_start, ev_total_stop, timing.total_ms);

    // Cleanup events
    cudaEventDestroy(ev_total_start);  cudaEventDestroy(ev_total_stop);
    cudaEventDestroy(ev_h2d_start);    cudaEventDestroy(ev_h2d_stop);
    cudaEventDestroy(ev_desc_start);   cudaEventDestroy(ev_desc_stop);
    cudaEventDestroy(ev_fr_start);     cudaEventDestroy(ev_fr_stop);
    cudaEventDestroy(ev_fa_start);     cudaEventDestroy(ev_fa_stop);
    cudaEventDestroy(ev_d2h_start);    cudaEventDestroy(ev_d2h_stop);

    // Cleanup GPU memory
    cudaFree(d_type);
    cudaFree(d_NN_r); cudaFree(d_NL_r);
    cudaFree(d_NN_a); cudaFree(d_NL_a);
    cudaFree(d_x12_r); cudaFree(d_y12_r); cudaFree(d_z12_r);
    cudaFree(d_x12_a); cudaFree(d_y12_a); cudaFree(d_z12_a);
    cudaFree(d_rc_r); cudaFree(d_rc_a);
    cudaFree(d_ann_c);
    cudaFree(d_w0); cudaFree(d_b0); cudaFree(d_w1); cudaFree(d_b1);
    cudaFree(d_qs);
    cudaFree(d_pot); cudaFree(d_Fp); cudaFree(d_sfxyz);
    cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz); cudaFree(d_vir);
}
