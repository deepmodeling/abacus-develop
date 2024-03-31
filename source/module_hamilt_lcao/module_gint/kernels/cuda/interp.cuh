#ifndef INTERP_CUH
#define INTERP_CUH

#include <cuda_runtime.h>

namespace GintKernel
{
    static __device__
    void interpolate(double distance, double delta_r_g, int it, double nwmax_g, int nr_max,
                     int *atom_nw, bool *atom_iw2_new, double *psi_u, double ylma[49], int *atom_iw2_ylm,
                     double *psir_ylm_left, int dist_tmp, int stride)
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
                phi = c1 * psi_u[iw_nr] + c2 * psi_u[iw_nr + 1] + c3 * psi_u[iw_nr + 2] + c4 * psi_u[iw_nr + 3];
            }
            psir_ylm_left[dist_tmp] = phi * ylma[atom_iw2_ylm[it_nw_iw]];
            dist_tmp += stride;
            iw_nr += 2 * nr_max;
            it_nw_iw++;
        }
    }
}

#endif