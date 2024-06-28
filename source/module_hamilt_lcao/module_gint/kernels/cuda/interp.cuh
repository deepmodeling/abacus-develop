#ifndef INTERP_CUH
#define INTERP_CUH

#include <cuda_runtime.h>

namespace GintKernel
{
static __device__ void interpolate(const double dist,
                                   const double delta_r_g,
                                   const int it,
                                   const double nwmax_g,
                                   const int nr_max,
                                   const int* __restrict__ atom_nw,
                                   const bool* __restrict__ atom_iw2_new,
                                   const double* __restrict__ psi_u,
                                   const double ylma[49],
                                   const int* __restrict__ atom_iw2_ylm,
                                   double* psir_ylm_left,
                                   int dist_tmp,
                                   const int stride)
{
    const double distance = dist / delta_r_g;

    const int ip = (int)(distance);
    const double dx = distance - ip;
    const double dx2 = dx * dx;
    const double dx3 = dx2 * dx;

    const double c3 = 3.0 * dx2 - 2.0 * dx3;
    const double c1 = 1.0 - c3;
    const double c2 = (dx - 2.0 * dx2 + dx3) * delta_r_g;
    const double c4 = (dx3 - dx2) * delta_r_g;

    double phi = 0.0;
    const int it_nw = it * nwmax_g;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    for (int iw = 0; iw < atom_nw[it]; ++iw)
    {
        if (atom_iw2_new[it_nw_iw])
        {
            phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1]
                  + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
        }
        psir_ylm_left[dist_tmp] = phi * ylma[atom_iw2_ylm[it_nw_iw]];
        dist_tmp += stride;
        iw_nr += 2 * nr_max;
        it_nw_iw++;
    }
}

static __device__ void interpolate_f(const double distance,
                                     const double delta_r_g,
                                     const int it,
                                     const int nwmax_g,
                                     const int nr_max,
                                     const int* __restrict__ atom_nw,
                                     const bool* __restrict__ atom_iw2_new,
                                     const double* __restrict__ psi_u,
                                     const int* __restrict__ atom_iw2_l,
                                     const int* __restrict__ atom_iw2_ylm,
                                     double* psi,
                                     int dist_tmp,
                                     const double ylma[49],
                                     const double vlbr3_value,
                                     double* dpsi,
                                     const double * __restrict__ dr,
                                     const double grly[49][3],
                                     double* d2psi)
{
    // Calculate normalized position for interpolation
    const double postion = distance / delta_r_g;
    // Extract integer part and fractional part of the position
    const double ip = static_cast<int>(postion);
    const double x0 = postion - ip;
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double x12 = x1 * x2 / 6;
    const double x03 = x0 * x3 / 2;
    // Temporary variables for interpolation
    double tmp = 0.0;
    double dtmp = 0.0;
    // Loop over non-zero elements in atom_nw array
    const int it_nw = it * nwmax_g;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    double dpsir[150][4]={0.0};
    int dist_tmp_clac=dist_tmp;
    for (int iw = 0; iw < atom_nw[it]; ++iw)
    {
        if (atom_iw2_new[it_nw_iw])
        {
            // Perform interpolation using cubic B-spline
            // basis functions
            tmp = x12 * (psi_u[iw_nr] * x3 + psi_u[iw_nr + 6] * x0)
                  + x03 * (psi_u[iw_nr + 2] * x2 - psi_u[iw_nr + 4] * x1);
            dtmp = x12 * (psi_u[iw_nr + 1] * x3 + psi_u[iw_nr + 7] * x0)
                   + x03 * (psi_u[iw_nr + 3] * x2 - psi_u[iw_nr + 5] * x1);
        }
        // Extract information from atom_iw2_* arrays
        const int ll = atom_iw2_l[it_nw_iw];

        const int idx_lm = atom_iw2_ylm[it_nw_iw];
        const double rl = pow(distance, ll);
        const double rl_r = 1.0 / rl;
        const double dist_r = 1 / distance;
        const int dist_tmp_force = dist_tmp_clac * 3;
        const int dist_tmp_stress = dist_tmp_clac * 6;
        // Compute right-hand side of the equation
        dpsir[iw][3] = tmp * ylma[idx_lm] * rl_r * vlbr3_value;
        // Compute derivatives with respect to spatial
        // coordinates
        const double tmpdphi_rly
            = (dtmp - tmp * ll * dist_r) * rl_r * ylma[idx_lm] * dist_r;
        const double tmprl = tmp * rl_r;
        const double dpsirx = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];
        const double dpsiry = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
        const double dpsirz = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];
        dpsir[iw][0] = dpsirx;
        dpsir[iw][1] = dpsiry;
        dpsir[iw][2] = dpsirz;

        // Update loop counters and indices
        dist_tmp_clac += 1;
        iw_nr += nr_max;
        iw_nr += nr_max;
        it_nw_iw++;
    }

    #pragma unroll
    int dist_tmp_trans = dist_tmp;
    for (int iw=0;iw<atom_nw[it];++iw)
    {
        const int dist_tmp_force = dist_tmp_trans * 3;
        const int dist_tmp_stress = dist_tmp_trans * 6;
        psi[dist_tmp_trans] = dpsir[iw][3];
        dpsi[dist_tmp_force] = dpsir[iw][0];
        dpsi[dist_tmp_force + 1] = dpsir[iw][1];
        dpsi[dist_tmp_force + 2] = dpsir[iw][2];

        d2psi[dist_tmp_stress] = dpsir[iw][0] * dr[0];
        d2psi[dist_tmp_stress + 1] = dpsir[iw][0] * dr[1];
        d2psi[dist_tmp_stress + 2] = dpsir[iw][0] * dr[2];
        d2psi[dist_tmp_stress + 3] = dpsir[iw][1] * dr[1];
        d2psi[dist_tmp_stress + 4] = dpsir[iw][1] * dr[2];
        d2psi[dist_tmp_stress + 5] = dpsir[iw][2] * dr[2];
        dist_tmp_trans += 1;
    }
}
} // namespace GintKernel

#endif