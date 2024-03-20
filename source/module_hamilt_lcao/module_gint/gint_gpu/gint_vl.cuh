#ifndef GINT_VL_CUH
#define GINT_VL_CUH

#include <cuda_runtime.h>
namespace lcaoCudaKernel{
__global__ void get_psi_and_vldr3(double *ylmcoef,
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
                                  double *psir_ylm_left,
                                  double *psir_ylm_right);

__global__ void psi_multiple(const int* m, int* n,
                                double  const * const * global_A_array,
                                double const * const * global_B_array,
                                double ** global_C_array);

#endif // GINT_VL_CUH
}