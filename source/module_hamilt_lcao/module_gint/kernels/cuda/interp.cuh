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
                                   int* atom_nw,
                                   bool* atom_iw2_new,
                                   double* psi_u,
                                   double ylma[49],
                                   int* atom_iw2_ylm,
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

static __device__ void interpolate_f(double& distance,
                                     double delta_r_g,
                                     int it,
                                     double nwmax_g,
                                     int nr_max,
                                     int* atom_nw,
                                     bool* atom_iw2_new,
                                     double* psi_u,
                                     int* atom_iw2_l,
                                     int* atom_iw2_ylm,
                                     double* psir_ylm_right,
                                     int& dist_tmp,
                                     double ylma[49],
                                     double vlbr3_value,
                                     double* dpsir_ylm_left_x,
                                     double dr[3],
                                     double grly[49][3],
                                     double* dpsir_ylm_left_y,
                                     double* dpsir_ylm_left_z,
                                     double* ddpsir_ylm_left_xx,
                                     double* ddpsir_ylm_left_xy,
                                     double* ddpsir_ylm_left_xz,
                                     double* ddpsir_ylm_left_yy,
                                     double* ddpsir_ylm_left_yz,
                                     double* ddpsir_ylm_left_zz)
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
        psir_ylm_right[dist_tmp] = tmp * ylma[idx_lm] / rl * vlbr3_value;
        // Compute derivatives with respect to spatial
        // coordinates
        const double tmpdphi_rly
            = (dtmp - tmp * ll / distance) / rl * ylma[idx_lm] / distance;
        const double tmprl = tmp / rl;
        dpsir_ylm_left_x[dist_tmp]
            = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];

        dpsir_ylm_left_y[dist_tmp]
            = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
        dpsir_ylm_left_z[dist_tmp]
            = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];

        ddpsir_ylm_left_xx[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[0];
        ddpsir_ylm_left_xy[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[1];
        ddpsir_ylm_left_xz[dist_tmp] = dpsir_ylm_left_x[dist_tmp] * dr[2];
        ddpsir_ylm_left_yy[dist_tmp] = dpsir_ylm_left_y[dist_tmp] * dr[1];
        ddpsir_ylm_left_yz[dist_tmp] = dpsir_ylm_left_y[dist_tmp] * dr[2];
        ddpsir_ylm_left_zz[dist_tmp] = dpsir_ylm_left_z[dist_tmp] * dr[2];

        // Update loop counters and indices
        dist_tmp += 1;
        iw_nr += nr_max;
        iw_nr += nr_max;
        it_nw_iw++;
    }
}
} // namespace GintKernel

#endif