#ifndef GINT_VL_CUH
#define GINT_VL_CUH

#include <cuda_runtime.h>
namespace GintKernel{
/*
 * @brief: get the value of the spherical harmonics
 * 
 * 
 * @note the left and right matrix elements of the grid point integral. 
 * We can understand the grid point integral of the local potential term 
 * as the following operation: 
 * H = psi * vlocal * psi * dr^3. 
 * Here, the matrix element of the left matrix is psi, and the matrix 
 * element of the right matrix is vlocal * psi * dr^3.
*/

__global__ void get_psi_and_vldr3(double *ylmcoef,
                                  double delta_r_g,
                                  int bxyz_g,
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

} // namespace GintKernel
#endif // GINT_VL_CUH