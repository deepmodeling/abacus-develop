#ifndef GINT_RHO_CUH
#define GINT_RHO_CUH

#include <cuda_runtime.h>
namespace lcaoCudaKernel{
__global__ void get_psi(double *ylmcoef,
                        double delta_r_g,
                        double bxyz_g,
                        double nwmax_g,
                        double *input_double,
                        int *input_int,
                        int *num_psir,
                        int psi_size_max,
                        int *ucell_atom_nwl,
                        bool *atom_iw2_new,
                        int *atom_iw2_ylm,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *psir_ylm_left);

__global__ void psir_dot(int * n,
                        double **x_array_g,
                        int incx,
                        double **y_array_g,
                        int incy,
                        double **results_g,
                        int batchcount);

#endif // GINT_RHO_CUH
}