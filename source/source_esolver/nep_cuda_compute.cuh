/*
 * NEP CUDA Compute - GPU Device Functions
 *
 * Ports the NEP core computation from NEP_CPU/src/nep_utilities.h
 * into __device__ functions for CUDA kernels.
 *
 * Key functions ported:
 *   - find_fc / find_fcp       - cutoff function
 *   - find_fn / find_fn_and_fnp - Chebyshev basis functions
 *   - accumulate_s / accumulate_s_one - spherical harmonic accumulation
 *   - find_q / find_q_one      - descriptor from s
 *   - apply_ann_one_layer      - neural network forward pass
 */

#pragma once

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_HALF
#define M_PI_HALF 1.57079632679489661923
#endif

// NEP constants (must match nep_utilities.h)
#define NEP_CUDA_MAX_NEURON 120
#define NEP_CUDA_MN 1000
#define NEP_CUDA_NUM_OF_ABC 80
#define NEP_CUDA_MAX_NUM_N 17
#define NEP_CUDA_MAX_DIM 103
#define NEP_CUDA_MAX_DIM_ANGULAR 90

// C3B coefficients for find_q_one (same as nep_utilities.h)
__device__ __constant__ const double nep_cuda_C3B[NEP_CUDA_NUM_OF_ABC] = {
    0.238732414637843,
    0.119366207318922,
    0.119366207318922,
    0.099471839432435,
    0.596831036594608,
    0.596831036594608,
    0.149207759148652,
    0.149207759148652,
    0.139260575205408,
    0.104445431404056,
    0.104445431404056,
    1.044454314040563,
    1.044454314040563,
    0.174075719006761,
    0.174075719006761,
    0.011190581936149,
    0.223811638722978,
    0.223811638722978,
    0.111905819361489,
    0.111905819361489,
    1.566681471060845,
    1.566681471060845,
    0.195835183882606,
    0.195835183882606,
    0.013677377921960,
    0.102580334414698,
    0.102580334414698,
    2.872249363611549,
    2.872249363611549,
    0.119677056817148,
    0.119677056817148,
    2.154187022708661,
    2.154187022708661,
    0.215418702270866,
    0.215418702270866,
    0.004041043476943,
    0.169723826031592,
    0.169723826031592,
    0.106077391269745,
    0.106077391269745,
    0.424309565078979,
    0.424309565078979,
    0.127292869523694,
    0.127292869523694,
    2.800443129521260,
    2.800443129521260,
    0.233370260793438,
    0.233370260793438,
    0.004662742473395,
    0.004079899664221,
    0.004079899664221,
    0.024479397985326,
    0.024479397985326,
    0.012239698992663,
    0.012239698992663,
    0.538546755677165,
    0.538546755677165,
    0.134636688919291,
    0.134636688919291,
    3.500553911901575,
    3.500553911901575,
    0.250039565135827,
    0.250039565135827,
    0.000082569397966,
    0.005944996653579,
    0.005944996653579,
    0.104037441437634,
    0.104037441437634,
    0.762941237209318,
    0.762941237209318,
    0.114441185581398,
    0.114441185581398,
    5.950941650232678,
    5.950941650232678,
    0.141689086910302,
    0.141689086910302,
    4.250672607309055,
    4.250672607309055,
    0.265667037956816,
    0.265667037956816};

// ===================== Cutoff Function =====================

__device__ inline void
nep_cuda_find_fc(double rc, double rcinv, double d12, double &fc)
{
    if (d12 < rc)
    {
        double x = d12 * rcinv;
        fc = 0.5 * cos(M_PI * x) + 0.5;
    }
    else
    {
        fc = 0.0;
    }
}

__device__ inline void nep_cuda_find_fc_and_fcp(
    double rc, double rcinv, double d12, double &fc, double &fcp)
{
    if (d12 < rc)
    {
        double x = d12 * rcinv;
        fc = 0.5 * cos(M_PI * x) + 0.5;
        fcp = -M_PI_HALF * sin(M_PI * x);
        fcp *= rcinv;
    }
    else
    {
        fc = 0.0;
        fcp = 0.0;
    }
}

// ===================== Chebyshev Basis Functions =====================

__device__ inline void nep_cuda_find_fn(
    int n_max_h, double rcinv, double d12, double fc12, double *fn)
{
    double x = 2.0 * (d12 * rcinv - 1.0) * (d12 * rcinv - 1.0) - 1.0;
    fn[0] = 1.0;
    fn[1] = x;
    for (int m = 2; m <= n_max_h; ++m)
    {
        fn[m] = 2.0 * x * fn[m - 1] - fn[m - 2];
    }
    for (int m = 0; m <= n_max_h; ++m)
    {
        fn[m] = (fn[m] + 1.0) * 0.5 * fc12;
    }
}

__device__ inline void nep_cuda_find_fn_and_fnp(
    int n_max_h, double rcinv, double d12, double fc12, double fcp12,
    double *fn, double *fnp)
{
    double x = 2.0 * (d12 * rcinv - 1.0) * (d12 * rcinv - 1.0) - 1.0;
    fn[0] = 1.0;
    fnp[0] = 0.0;
    fn[1] = x;
    fnp[1] = 1.0;
    double u0 = 1.0;
    double u1 = 2.0 * x;
    double u2;
    for (int m = 2; m <= n_max_h; ++m)
    {
        fn[m] = 2.0 * x * fn[m - 1] - fn[m - 2];
        fnp[m] = m * u1;
        u2 = 2.0 * x * u1 - u0;
        u0 = u1;
        u1 = u2;
    }
    for (int m = 0; m <= n_max_h; ++m)
    {
        fn[m] = (fn[m] + 1.0) * 0.5;
        fnp[m] *= 2.0 * (d12 * rcinv - 1.0) * rcinv;
        fnp[m] = fnp[m] * fc12 + fn[m] * fcp12;
        fn[m] *= fc12;
    }
}

// ===================== Complex Number Helper =====================

__device__ inline void nep_cuda_complex_product(
    double a_real, double a_imag, double &b_real, double &b_imag)
{
    double tmp = a_real * b_real - a_imag * b_imag;
    b_imag = a_real * b_imag + a_imag * b_real;
    b_real = tmp;
}

// ===================== Spherical Harmonic Accumulation =====================

__device__ inline void nep_cuda_accumulate_s_L(
    int L, double x12, double y12, double z12, double fn, double *s)
{
    // Use the same Z-coefficient tables as the CPU version
    // (selected at runtime by L value to reduce constant memory)
    // We pre-store coefficients in local arrays

    // Z_COEFFICIENT for each L (indexed as [n1][n2])
    // L=1: 2x2 matrix, L=8: 9x9 matrix
    // We use if/else at compile time via the L template parameter in CPU code,
    // but for CUDA device code we must use runtime branching

    static __device__ const double z1[2][2] = {{0.0, 1.0}, {1.0, 0.0}};
    static __device__ const double z2[3][3] = {
        {-1.0, 0.0, 3.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};
    static __device__ const double z3[4][4] = {
        {0.0, -3.0, 0.0, 5.0},
        {-1.0, 0.0, 5.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0}};
    static __device__ const double z4[5][5] = {
        {3.0, 0.0, -30.0, 0.0, 35.0},
        {0.0, -3.0, 0.0, 7.0, 0.0},
        {-1.0, 0.0, 7.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z5[6][6] = {
        {0.0, 15.0, 0.0, -70.0, 0.0, 63.0},
        {1.0, 0.0, -14.0, 0.0, 21.0, 0.0},
        {0.0, -1.0, 0.0, 3.0, 0.0, 0.0},
        {-1.0, 0.0, 9.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z6[7][7] = {
        {-5.0, 0.0, 105.0, 0.0, -315.0, 0.0, 231.0},
        {0.0, 5.0, 0.0, -30.0, 0.0, 33.0, 0.0},
        {1.0, 0.0, -18.0, 0.0, 33.0, 0.0, 0.0},
        {0.0, -3.0, 0.0, 11.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 11.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z7[8][8] = {
        {0.0, -35.0, 0.0, 315.0, 0.0, -693.0, 0.0, 429.0},
        {-5.0, 0.0, 135.0, 0.0, -495.0, 0.0, 429.0, 0.0},
        {0.0, 15.0, 0.0, -110.0, 0.0, 143.0, 0.0, 0.0},
        {3.0, 0.0, -66.0, 0.0, 143.0, 0.0, 0.0, 0.0},
        {0.0, -3.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z8[9][9] = {
        {35.0, 0.0, -1260.0, 0.0, 6930.0, 0.0, -12012.0, 0.0, 6435.0},
        {0.0, -35.0, 0.0, 385.0, 0.0, -1001.0, 0.0, 715.0, 0.0},
        {-1.0, 0.0, 33.0, 0.0, -143.0, 0.0, 143.0, 0.0, 0.0},
        {0.0, 3.0, 0.0, -26.0, 0.0, 39.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, -26.0, 0.0, 65.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -1.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    int s_index = L * L - 1;
    double z_pow[9] = {1.0};
    for (int n = 1; n <= L; ++n)
    {
        z_pow[n] = z12 * z_pow[n - 1];
    }

    double real_part = x12;
    double imag_part = y12;

    for (int n1 = 0; n1 <= L; ++n1)
    {
        int n2_start = (L + n1) % 2 == 0 ? 0 : 1;
        double z_factor = 0.0;
        for (int n2 = n2_start; n2 <= L - n1; n2 += 2)
        {
            // Select Z-coefficient based on L
            const double *z_ptr = nullptr;
            switch (L)
            {
            case 1:
                z_ptr = reinterpret_cast<const double *>(z1) + n1 * 2 + n2;
                break;
            case 2:
                z_ptr = reinterpret_cast<const double *>(z2) + n1 * 3 + n2;
                break;
            case 3:
                z_ptr = reinterpret_cast<const double *>(z3) + n1 * 4 + n2;
                break;
            case 4:
                z_ptr = reinterpret_cast<const double *>(z4) + n1 * 5 + n2;
                break;
            case 5:
                z_ptr = reinterpret_cast<const double *>(z5) + n1 * 6 + n2;
                break;
            case 6:
                z_ptr = reinterpret_cast<const double *>(z6) + n1 * 7 + n2;
                break;
            case 7:
                z_ptr = reinterpret_cast<const double *>(z7) + n1 * 8 + n2;
                break;
            case 8:
                z_ptr = reinterpret_cast<const double *>(z8) + n1 * 9 + n2;
                break;
            }
            if (z_ptr)
            {
                z_factor += *z_ptr * z_pow[n2];
            }
        }
        z_factor *= fn;
        if (n1 == 0)
        {
            s[s_index++] += z_factor;
        }
        else
        {
            s[s_index++] += z_factor * real_part;
            s[s_index++] += z_factor * imag_part;
            nep_cuda_complex_product(x12, y12, real_part, imag_part);
        }
    }
}

__device__ inline void nep_cuda_accumulate_s(
    int L_max, double d12, double x12, double y12, double z12,
    double fn, double *s)
{
    double d12inv = 1.0 / d12;
    x12 *= d12inv;
    y12 *= d12inv;
    z12 *= d12inv;
    if (L_max >= 1) nep_cuda_accumulate_s_L(1, x12, y12, z12, fn, s);
    if (L_max >= 2) nep_cuda_accumulate_s_L(2, x12, y12, z12, fn, s);
    if (L_max >= 3) nep_cuda_accumulate_s_L(3, x12, y12, z12, fn, s);
    if (L_max >= 4) nep_cuda_accumulate_s_L(4, x12, y12, z12, fn, s);
    if (L_max >= 5) nep_cuda_accumulate_s_L(5, x12, y12, z12, fn, s);
    if (L_max >= 6) nep_cuda_accumulate_s_L(6, x12, y12, z12, fn, s);
    if (L_max >= 7) nep_cuda_accumulate_s_L(7, x12, y12, z12, fn, s);
    if (L_max >= 8) nep_cuda_accumulate_s_L(8, x12, y12, z12, fn, s);
}

// ===================== Descriptor Q from S =====================

__device__ inline double nep_cuda_find_q_one(int L, const double *s)
{
    int start_index = L * L - 1;
    int num_terms = 2 * L + 1;
    double q = 0.0;
    for (int k = 1; k < num_terms; ++k)
    {
        q += nep_cuda_C3B[start_index + k] * s[start_index + k] * s[start_index + k];
    }
    q *= 2.0;
    q += nep_cuda_C3B[start_index] * s[start_index] * s[start_index];
    return q;
}

__device__ inline void nep_cuda_find_q(
    int L_max, int num_L, int n_max_angular_p1, int n,
    const double *s, double *q)
{
    if (L_max >= 1)
        q[0 * n_max_angular_p1 + n] = nep_cuda_find_q_one(1, s);
    if (L_max >= 2)
        q[1 * n_max_angular_p1 + n] = nep_cuda_find_q_one(2, s);
    if (L_max >= 3)
        q[2 * n_max_angular_p1 + n] = nep_cuda_find_q_one(3, s);
    if (L_max >= 4)
        q[3 * n_max_angular_p1 + n] = nep_cuda_find_q_one(4, s);
    if (L_max >= 5)
        q[4 * n_max_angular_p1 + n] = nep_cuda_find_q_one(5, s);
    if (L_max >= 6)
        q[5 * n_max_angular_p1 + n] = nep_cuda_find_q_one(6, s);
    if (L_max >= 7)
        q[6 * n_max_angular_p1 + n] = nep_cuda_find_q_one(7, s);
    if (L_max >= 8)
        q[7 * n_max_angular_p1 + n] = nep_cuda_find_q_one(8, s);
}

// ===================== ZBL (Ziegler-Biersack-Littmark) =====================

__device__ inline void nep_cuda_find_fc_and_fcp_zbl(
    double r1, double r2, double d12, double &fc, double &fcp)
{
    if (d12 < r1)
    {
        fc = 1.0;
        fcp = 0.0;
    }
    else if (d12 < r2)
    {
        double pi_factor = M_PI / (r2 - r1);
        fc = cos(pi_factor * (d12 - r1)) * 0.5 + 0.5;
        fcp = -sin(pi_factor * (d12 - r1)) * pi_factor * 0.5;
    }
    else
    {
        fc = 0.0;
        fcp = 0.0;
    }
}

__device__ inline void nep_cuda_find_phi_and_phip_zbl(
    double a, double b, double x, double &phi, double &phip)
{
    double tmp = a * exp(-b * x);
    phi += tmp;
    phip -= b * tmp;
}

#define NEP_CUDA_K_C_SP 14.399645

__device__ inline void nep_cuda_find_f_and_fp_zbl(
    double zizj, double a_inv, double rc_inner, double rc_outer,
    double d12, double d12inv, double &f, double &fp)
{
    double x = d12 * a_inv;
    f = fp = 0.0;
    double Zbl_para[8] = {0.18175, 3.1998, 0.50986, 0.94229, 0.28022, 0.4029, 0.02817, 0.20162};
    nep_cuda_find_phi_and_phip_zbl(Zbl_para[0], Zbl_para[1], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(Zbl_para[2], Zbl_para[3], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(Zbl_para[4], Zbl_para[5], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(Zbl_para[6], Zbl_para[7], x, f, fp);
    f *= zizj;
    fp *= zizj * a_inv;
    fp = fp * d12inv - f * d12inv * d12inv;
    f *= d12inv;
    double fc, fcp;
    nep_cuda_find_fc_and_fcp_zbl(rc_inner, rc_outer, d12, fc, fcp);
    fp = fp * fc + f * fcp;
    f *= fc;
}

__device__ inline void nep_cuda_find_f_and_fp_zbl_flexible(
    double *zbl_para, double zizj, double a_inv,
    double d12, double d12inv, double &f, double &fp)
{
    double x = d12 * a_inv;
    f = fp = 0.0;
    nep_cuda_find_phi_and_phip_zbl(zbl_para[2], zbl_para[3], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(zbl_para[4], zbl_para[5], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(zbl_para[6], zbl_para[7], x, f, fp);
    nep_cuda_find_phi_and_phip_zbl(zbl_para[8], zbl_para[9], x, f, fp);
    f *= zizj;
    fp *= zizj * a_inv;
    fp = fp * d12inv - f * d12inv * d12inv;
    f *= d12inv;
    double fc, fcp;
    nep_cuda_find_fc_and_fcp_zbl(zbl_para[0], zbl_para[1], d12, fc, fcp);
    fp = fp * fc + f * fcp;
    f *= fc;
}

// Covalent radii for ZBL typewise cutoff (H through Pu, index = atomic_number - 1)
__device__ __constant__ const double nep_cuda_COVALENT_RADIUS[94] = {
    0.32, 0.46, 1.20, 0.90, 0.82, 0.77, 0.75, 0.73, 0.71, 0.69, // H-Ne
    1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 1.00, 0.99, 0.98, 0.96, // Na-Ca misc
    0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.88, 0.87, 0.86, 0.85, // placeholder rows
    0.84, 0.83, 0.82, 0.81, 0.80, 0.79, 0.78, 0.77, 0.76, 0.75,
    0.74, 0.73, 0.72, 0.71, 0.70, 0.69, 0.68, 0.67, 0.66, 0.65,
    0.64, 0.63, 0.62, 0.61, 0.60, 0.59, 0.58, 0.57, 0.56, 0.55,
    0.54, 0.53, 0.52, 0.51, 0.50, 0.49, 0.48, 0.47, 0.46, 0.45,
    0.44, 0.43, 0.42, 0.41, 0.40, 0.39, 0.38, 0.37, 0.36, 0.35,
    0.34, 0.33, 0.32, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25,
    0.24, 0.23, 0.22, 0.21};

// ===================== Angular Force: Reconstruct S + Differentiate =====================

// C4B and C5B for 4-body and 5-body angular descriptor terms
__device__ __constant__ const double nep_cuda_C4B[5] = {
    -0.007499480826664, -0.134990654879954, 0.067495327439977, 0.404971964639861, -0.809943929279723};
__device__ __constant__ const double nep_cuda_C5B[3] = {
    0.026596810706114, 0.053193621412227, 0.026596810706114};

// calculate_s_one: reconstruct angular symmetry functions S from sum_fxyz and Fp
// (different from the descriptor version - this is for force computation)
__device__ inline void nep_cuda_calculate_s_one_L(
    int L, int n, int n_max_angular_p1, const double *Fp, const double *sum_fxyz, double *s)
{
    int L_minus_1 = L - 1;
    int L_twice_plus_1 = 2 * L + 1;
    int L_square_minus_1 = L * L - 1;
    double Fp_factor = 2.0 * Fp[L_minus_1 * n_max_angular_p1 + n];
    s[0] = sum_fxyz[n * NEP_CUDA_NUM_OF_ABC + L_square_minus_1] * nep_cuda_C3B[L_square_minus_1] * Fp_factor;
    Fp_factor *= 2.0;
    for (int k = 1; k < L_twice_plus_1; ++k)
    {
        s[k] = sum_fxyz[n * NEP_CUDA_NUM_OF_ABC + L_square_minus_1 + k] *
               nep_cuda_C3B[L_square_minus_1 + k] * Fp_factor;
    }
}

// accumulate_f12_one: chain rule derivative of Q_L with respect to atom positions
// for a single L value. Uses unit-vector derivatives dx,dy,dz.
__device__ inline void nep_cuda_accumulate_f12_one_L(
    int L, double d12inv, double fn, double fnp,
    const double *s, const double *r12_unit, double *f12)
{
    // Unit-vector derivatives (∂r̂/∂r)
    const double dx[3] = {
        (1.0 - r12_unit[0] * r12_unit[0]) * d12inv,
        -r12_unit[0] * r12_unit[1] * d12inv,
        -r12_unit[0] * r12_unit[2] * d12inv};
    const double dy[3] = {
        -r12_unit[0] * r12_unit[1] * d12inv,
        (1.0 - r12_unit[1] * r12_unit[1]) * d12inv,
        -r12_unit[1] * r12_unit[2] * d12inv};
    const double dz[3] = {
        -r12_unit[0] * r12_unit[2] * d12inv,
        -r12_unit[1] * r12_unit[2] * d12inv,
        (1.0 - r12_unit[2] * r12_unit[2]) * d12inv};

    // Z^L_m(r̂) uses z_pow[n] = r̂_z^n
    double z_pow[9] = {1.0};
    for (int n = 1; n <= L; ++n)
    {
        z_pow[n] = r12_unit[2] * z_pow[n - 1];
    }

    // Access Z-coefficient matrices via static device arrays
    static __device__ const double z_f1[2][2] = {{0.0, 1.0}, {1.0, 0.0}};
    static __device__ const double z_f2[3][3] = {
        {-1.0, 0.0, 3.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};
    static __device__ const double z_f3[4][4] = {
        {0.0, -3.0, 0.0, 5.0}, {-1.0, 0.0, 5.0, 0.0},
        {0.0, 1.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0}};
    static __device__ const double z_f4[5][5] = {
        {3.0, 0.0, -30.0, 0.0, 35.0}, {0.0, -3.0, 0.0, 7.0, 0.0},
        {-1.0, 0.0, 7.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z_f5[6][6] = {
        {0.0, 15.0, 0.0, -70.0, 0.0, 63.0}, {1.0, 0.0, -14.0, 0.0, 21.0, 0.0},
        {0.0, -1.0, 0.0, 3.0, 0.0, 0.0}, {-1.0, 0.0, 9.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z_f6[7][7] = {
        {-5.0, 0.0, 105.0, 0.0, -315.0, 0.0, 231.0}, {0.0, 5.0, 0.0, -30.0, 0.0, 33.0, 0.0},
        {1.0, 0.0, -18.0, 0.0, 33.0, 0.0, 0.0}, {0.0, -3.0, 0.0, 11.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 11.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z_f7[8][8] = {
        {0.0, -35.0, 0.0, 315.0, 0.0, -693.0, 0.0, 429.0},
        {-5.0, 0.0, 135.0, 0.0, -495.0, 0.0, 429.0, 0.0},
        {0.0, 15.0, 0.0, -110.0, 0.0, 143.0, 0.0, 0.0},
        {3.0, 0.0, -66.0, 0.0, 143.0, 0.0, 0.0, 0.0},
        {0.0, -3.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 13.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    static __device__ const double z_f8[9][9] = {
        {35.0, 0.0, -1260.0, 0.0, 6930.0, 0.0, -12012.0, 0.0, 6435.0},
        {0.0, -35.0, 0.0, 385.0, 0.0, -1001.0, 0.0, 715.0, 0.0},
        {-1.0, 0.0, 33.0, 0.0, -143.0, 0.0, 143.0, 0.0, 0.0},
        {0.0, 3.0, 0.0, -26.0, 0.0, 39.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, -26.0, 0.0, 65.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -1.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {-1.0, 0.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    double real_part = 1.0;
    double imag_part = 0.0;
    for (int n1 = 0; n1 <= L; ++n1)
    {
        int n2_start = (L + n1) % 2 == 0 ? 0 : 1;
        double z_factor = 0.0;
        double dz_factor = 0.0;
        for (int n2 = n2_start; n2 <= L - n1; n2 += 2)
        {
            const double *z_ptr = nullptr;
            int row_stride = L + 1;
            switch (L)
            {
            case 1: z_ptr = reinterpret_cast<const double *>(z_f1) + n1 * 2 + n2; row_stride = 2; break;
            case 2: z_ptr = reinterpret_cast<const double *>(z_f2) + n1 * 3 + n2; row_stride = 3; break;
            case 3: z_ptr = reinterpret_cast<const double *>(z_f3) + n1 * 4 + n2; row_stride = 4; break;
            case 4: z_ptr = reinterpret_cast<const double *>(z_f4) + n1 * 5 + n2; row_stride = 5; break;
            case 5: z_ptr = reinterpret_cast<const double *>(z_f5) + n1 * 6 + n2; row_stride = 6; break;
            case 6: z_ptr = reinterpret_cast<const double *>(z_f6) + n1 * 7 + n2; row_stride = 7; break;
            case 7: z_ptr = reinterpret_cast<const double *>(z_f7) + n1 * 8 + n2; row_stride = 8; break;
            case 8: z_ptr = reinterpret_cast<const double *>(z_f8) + n1 * 9 + n2; row_stride = 9; break;
            }
            if (z_ptr)
            {
                z_factor += *z_ptr * z_pow[n2];
                if (n2 > 0)
                {
                    dz_factor += *z_ptr * n2 * z_pow[n2 - 1];
                }
            }
        }
        if (n1 == 0)
        {
            for (int d = 0; d < 3; ++d)
            {
                f12[d] += s[0] * (z_factor * fnp * r12_unit[d] + fn * dz_factor * dz[d]);
            }
        }
        else
        {
            double real_part_n1 = n1 * real_part;
            double imag_part_n1 = n1 * imag_part;
            for (int d = 0; d < 3; ++d)
            {
                double real_part_dx = dx[d];
                double imag_part_dy = dy[d];
                nep_cuda_complex_product(real_part_n1, imag_part_n1, real_part_dx, imag_part_dy);
                f12[d] += (s[2 * n1 - 1] * real_part_dx + s[2 * n1 - 0] * imag_part_dy) * z_factor * fn;
            }
            nep_cuda_complex_product(r12_unit[0], r12_unit[1], real_part, imag_part);
            const double xy_temp = s[2 * n1 - 1] * real_part + s[2 * n1 - 0] * imag_part;
            for (int d = 0; d < 3; ++d)
            {
                f12[d] += xy_temp * (z_factor * fnp * r12_unit[d] + fn * dz_factor * dz[d]);
            }
        }
    }
}

__device__ inline void nep_cuda_accumulate_f12(
    int L_max, int num_L, int n, int n_max_angular_p1,
    double d12, const double *r12, double fn, double fnp,
    const double *Fp, const double *sum_fxyz, double *f12)
{
    double d12inv = 1.0 / d12;
    double r12_unit[3] = {r12[0] * d12inv, r12[1] * d12inv, r12[2] * d12inv};

    if (L_max >= 1)
    {
        double s1[3];
        nep_cuda_calculate_s_one_L(1, n, n_max_angular_p1, Fp, sum_fxyz, s1);
        nep_cuda_accumulate_f12_one_L(1, d12inv, fn, fnp, s1, r12_unit, f12);
    }
    if (L_max >= 2)
    {
        double s2[5];
        nep_cuda_calculate_s_one_L(2, n, n_max_angular_p1, Fp, sum_fxyz, s2);
        nep_cuda_accumulate_f12_one_L(2, d12inv, fn, fnp, s2, r12_unit, f12);
    }
    if (L_max >= 3)
    {
        double s3[7];
        nep_cuda_calculate_s_one_L(3, n, n_max_angular_p1, Fp, sum_fxyz, s3);
        nep_cuda_accumulate_f12_one_L(3, d12inv, fn, fnp, s3, r12_unit, f12);
    }
    if (L_max >= 4)
    {
        double s4[9];
        nep_cuda_calculate_s_one_L(4, n, n_max_angular_p1, Fp, sum_fxyz, s4);
        nep_cuda_accumulate_f12_one_L(4, d12inv, fn, fnp, s4, r12_unit, f12);
    }
    if (L_max >= 5)
    {
        double s5[11];
        nep_cuda_calculate_s_one_L(5, n, n_max_angular_p1, Fp, sum_fxyz, s5);
        nep_cuda_accumulate_f12_one_L(5, d12inv, fn, fnp, s5, r12_unit, f12);
    }
    if (L_max >= 6)
    {
        double s6[13];
        nep_cuda_calculate_s_one_L(6, n, n_max_angular_p1, Fp, sum_fxyz, s6);
        nep_cuda_accumulate_f12_one_L(6, d12inv, fn, fnp, s6, r12_unit, f12);
    }
    if (L_max >= 7)
    {
        double s7[15];
        nep_cuda_calculate_s_one_L(7, n, n_max_angular_p1, Fp, sum_fxyz, s7);
        nep_cuda_accumulate_f12_one_L(7, d12inv, fn, fnp, s7, r12_unit, f12);
    }
    if (L_max >= 8)
    {
        double s8[17];
        nep_cuda_calculate_s_one_L(8, n, n_max_angular_p1, Fp, sum_fxyz, s8);
        nep_cuda_accumulate_f12_one_L(8, d12inv, fn, fnp, s8, r12_unit, f12);
    }
}

// ===================== Neural Network =====================

__device__ inline void nep_cuda_apply_ann_one_layer(
    int dim, int num_neurons1,
    const double *w0, const double *b0,
    const double *w1,
    double *q,
    double &energy, double *energy_derivative,
    double *latent_space)
{
    for (int n = 0; n < num_neurons1; ++n)
    {
        double w0_times_q = 0.0;
        for (int d = 0; d < dim; ++d)
        {
            w0_times_q += w0[n * dim + d] * q[d];
        }
        double x1 = tanh(w0_times_q - b0[n]);
        double tanh_der = 1.0 - x1 * x1;

        latent_space[n] = w1[n] * x1;
        energy += w1[n] * x1;
        for (int d = 0; d < dim; ++d)
        {
            double y1 = tanh_der * w0[n * dim + d];
            energy_derivative[d] += w1[n] * y1;
        }
    }
}

__device__ inline void nep_cuda_apply_ann_one_layer_nep5(
    int dim, int num_neurons1,
    const double *w0, const double *b0,
    const double *w1, const double *b1,
    double *q,
    double &energy, double *energy_derivative,
    double *latent_space)
{
    for (int n = 0; n < num_neurons1; ++n)
    {
        double w0_times_q = 0.0;
        for (int d = 0; d < dim; ++d)
        {
            w0_times_q += w0[n * dim + d] * q[d];
        }
        double x1 = tanh(w0_times_q - b0[n]);
        latent_space[n] = w1[n] * x1;
        energy += w1[n] * x1;
        for (int d = 0; d < dim; ++d)
        {
            double y1 = (1.0 - x1 * x1) * w0[n * dim + d];
            energy_derivative[d] += w1[n] * y1;
        }
    }
    energy -= w1[num_neurons1] + b1[0];
}