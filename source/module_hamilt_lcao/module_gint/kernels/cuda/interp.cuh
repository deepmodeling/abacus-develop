#ifndef INTERP_CUH
#define INTERP_CUH

#include <cuda_runtime.h>

namespace GintKernel
{
static __device__ void interpolate(double distance,
                                   double delta_r_g,
                                   int it,
                                   double nwmax_g,
                                   int nr_max,
                                   const int* const atom_nw,
                                   const bool* const atom_iw2_new,
                                   const double* const psi_u,
                                   const double ylma[49],
                                   const int* const atom_iw2_ylm,
                                   double* psir_ylm_left,
                                   int dist_tmp,
                                   int stride)
{
    distance /= delta_r_g;

    int ip = (int)(distance);
    double dx = distance - ip;
    double dx2 = dx * dx;
    double dx3 = dx2 * dx;

    double c3 = 3.0 * dx2 - 2.0 * dx3;
    double c1 = 1.0 - c3;
    double c2 = (dx - 2.0 * dx2 + dx3) * delta_r_g;
    double c4 = (dx3 - dx2) * delta_r_g;

    double phi = 0.0;
    int it_nw = it * nwmax_g;
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

static __device__ void interpolate_f(double distance,
                                     double delta_r_g,
                                     int it,
                                     double nwmax_g,
                                     int nr_max,
                                     const int* const atom_nw,
                                     const bool* const atom_iw2_new,
                                     const double* const psi_u,
                                     const int* const atom_iw2_l,
                                     const int* const atom_iw2_ylm,
                                     double* psir_r,
                                     int dist_tmp,
                                     const double ylma[49],
                                     double vlbr3_value,
                                     double* psir_lx,
                                     const double dr[3],
                                     const double grly[49][3],
                                     double* psir_ly,
                                     double* psir_lz,
                                     double* psir_lxx,
                                     double* psir_lxy,
                                     double* psir_lxz,
                                     double* psir_lyy,
                                     double* psir_lyz,
                                     double* psir_lzz)
{
    // Calculate normalized position for interpolation
    distance = sqrt(distance);
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
    double tmp, dtmp;
    // Loop over non-zero elements in atom_nw array
    int it_nw = it * nwmax_g;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
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

        // Compute right-hand side of the equation
        psir_r[dist_tmp] = tmp * ylma[idx_lm] / rl * vlbr3_value;
        // Compute derivatives with respect to spatial
        // coordinates
        const double tmpdphi_rly
            = (dtmp - tmp * ll / distance) / rl * ylma[idx_lm] / distance;
        const double tmprl = tmp / rl;
        psir_lx[dist_tmp]
            = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];

        psir_ly[dist_tmp]
            = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
        psir_lz[dist_tmp]
            = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];

        psir_lxx[dist_tmp] = psir_lx[dist_tmp] * dr[0];
        psir_lxy[dist_tmp] = psir_lx[dist_tmp] * dr[1];
        psir_lxz[dist_tmp] = psir_lx[dist_tmp] * dr[2];
        psir_lyy[dist_tmp] = psir_ly[dist_tmp] * dr[1];
        psir_lyz[dist_tmp] = psir_ly[dist_tmp] * dr[2];
        psir_lzz[dist_tmp] = psir_lz[dist_tmp] * dr[2];

        // Update loop counters and indices
        dist_tmp += 1;
        iw_nr += nr_max;
        iw_nr += nr_max;
        it_nw_iw++;
    }
}
} // namespace GintKernel

#endif