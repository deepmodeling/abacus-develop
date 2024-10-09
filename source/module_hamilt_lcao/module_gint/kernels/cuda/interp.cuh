#ifndef INTERP_CUH
#define INTERP_CUH

#include <cuda_runtime.h>

namespace GintKernel
{
// if exponent is an integer between 0 and 5 (the most common cases in gint),
// pow_int is much faster than std::pow
static __device__ double pow_int(double base, int exp)
{
    switch (exp)
    {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    default:
        double result = pow(base, exp);
        return result;      
    }
}

static __device__ void interp_rho(const double dist,
                                  const double delta_r,
                                  const int atype,
                                  const double nwmax,
                                  const int nr_max,
                                  const int* __restrict__ atom_nw,
                                  const bool* __restrict__ atom_iw2_new,
                                  const double* __restrict__ psi_u,
                                  const double ylma[49],
                                  const int* __restrict__ atom_iw2_ylm,
                                  double* psi,
                                  int psi_idx)
{
    const double distance = dist / delta_r;

    const int ip = (int)(distance);
    const double dx = distance - ip;
    const double dx2 = dx * dx;
    const double dx3 = dx2 * dx;

    const double c3 = 3.0 * dx2 - 2.0 * dx3;
    const double c1 = 1.0 - c3;
    const double c2 = (dx - 2.0 * dx2 + dx3) * delta_r;
    const double c4 = (dx3 - dx2) * delta_r;

    double phi = 0.0;
    const int it_nw = atype * nwmax;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    for (int iw = 0; iw < atom_nw[atype]; ++iw)
    {
        if (atom_iw2_new[it_nw_iw])
        {
            phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1]
                  + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
        }
        psi[psi_idx] = phi * ylma[atom_iw2_ylm[it_nw_iw]];
        psi_idx += 1;
        iw_nr += 2 * nr_max;
        it_nw_iw++;
    }
}

static __device__ void interp_vl(const double dist,
                                 const double delta_r,
                                 const int atype,
                                 const double nwmax,
                                 const int bxyz,
                                 const int nr_max,
                                 const int* __restrict__ atom_nw,
                                 const bool* __restrict__ atom_iw2_new,
                                 const double* __restrict__ psi_u,
                                 const double ylma[49],
                                 const int* __restrict__ atom_iw2_ylm,
                                 const double vldr3_value,
                                 double* psi,
                                 double* psi_vldr3,
                                 int psi_idx)
{
    const double distance = dist / delta_r;

    const int ip = (int)(distance);
    const double dx = distance - ip;
    const double dx2 = dx * dx;
    const double dx3 = dx2 * dx;

    const double c3 = 3.0 * dx2 - 2.0 * dx3;
    const double c1 = 1.0 - c3;
    const double c2 = (dx - 2.0 * dx2 + dx3) * delta_r;
    const double c4 = (dx3 - dx2) * delta_r;

    double phi = 0.0;
    const int it_nw = atype * nwmax;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    for (int iw = 0; iw < atom_nw[atype]; ++iw)
    {
        if (atom_iw2_new[it_nw_iw])
        {
            phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1]
                  + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
        }
        psi[psi_idx] = phi * ylma[atom_iw2_ylm[it_nw_iw]];
        psi_vldr3[psi_idx] = psi[psi_idx] * vldr3_value;
        psi_idx += bxyz;
        iw_nr += 2 * nr_max;
        it_nw_iw++;
    }
}

static __device__ void interp_f(const double dist,
                                const double delta_r,
                                const int atype,
                                const double nwmax,
                                const int nr_max,
                                const int* __restrict__ atom_nw,
                                const bool* __restrict__ atom_iw2_new,
                                const double* __restrict__ psi_u,
                                const double ylma[49],
                                const int* __restrict__ atom_iw2_l,
                                const int* __restrict__ atom_iw2_ylm,
                                const double vldr3_value,
                                const double * __restrict__ dr,
                                const double grly[49][3],
                                int psi_idx,
                                double* psi,
                                double* dpsi,
                                double* d2psi)
{
    // Calculate normalized position for interpolation
    const double postion = dist / delta_r;
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
    const int it_nw = atype * nwmax;
    int iw_nr = (it_nw * nr_max + ip) * 2;
    int it_nw_iw = it_nw;
    for (int iw = 0; iw < atom_nw[atype]; ++iw)
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
        const double rl = pow_int(dist, ll);
        const double rl_r = 1.0 / rl;
        const double dist_r = 1 / dist;
        const int dpsi_idx = psi_idx * 3;
        const int d2psi_idx = psi_idx * 6;
        // Compute derivatives with respect to spatial
        // coordinates
        const double tmpdphi_rly
            = (dtmp - tmp * ll * dist_r) * rl_r * ylma[idx_lm] * dist_r;
        const double tmprl = tmp * rl_r;
        const double dpsirx = tmpdphi_rly * dr[0] + tmprl * grly[idx_lm][0];
        const double dpsiry = tmpdphi_rly * dr[1] + tmprl * grly[idx_lm][1];
        const double dpsirz = tmpdphi_rly * dr[2] + tmprl * grly[idx_lm][2];

        psi[psi_idx] = tmprl * ylma[idx_lm] * vldr3_value;
        dpsi[dpsi_idx] = dpsirx;
        dpsi[dpsi_idx + 1] = dpsiry;
        dpsi[dpsi_idx + 2] = dpsirz;
        d2psi[d2psi_idx] = dpsirx * dr[0];
        d2psi[d2psi_idx + 1] = dpsirx * dr[1];
        d2psi[d2psi_idx + 2] = dpsirx * dr[2];
        d2psi[d2psi_idx + 3] = dpsiry * dr[1];
        d2psi[d2psi_idx + 4] = dpsiry * dr[2];
        d2psi[d2psi_idx + 5] = dpsirz * dr[2];
        // Update loop counters and indices
        psi_idx += 1;
        iw_nr += 2 * nr_max;
        it_nw_iw++;
    }
}
} // namespace GintKernel

#endif