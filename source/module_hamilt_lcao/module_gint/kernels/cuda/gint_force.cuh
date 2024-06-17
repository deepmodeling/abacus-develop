#ifndef GINT_FORCE_CUH
#define GINT_FORCE_CUH

#include <cuda_runtime.h>
namespace GintKernel
{

__global__ void get_psi_force(double* ylmcoef,
                              double delta_r_g,
                              int bxyz_g,
                              double nwmax_g,
                              double* __restrict__ psi_input_double,
                              int* __restrict__ psi_input_int,
                              int* __restrict__ atom_num_per_bcell,
                              int* __restrict__ start_idx_per_bcell,
                              const int* __restrict__ ucell_atom_nwl,
                              const bool* __restrict__ atom_iw2_new,
                              const int* __restrict__ atom_iw2_ylm,
                              const int* __restrict__ atom_iw2_l,
                              const int* __restrict__ atom_nw,
                              int nr_max,
                              const double* __restrict__ psi_u,
                              double* psi,
                              double* dpsi_dx,
                              double* dpsi_dy,
                              double* dpsi_dz,
                              double* d2psi_dxx,
                              double* d2psi_dxy,
                              double* d2psi_dxz,
                              double* d2psi_dyy,
                              double* d2psi_dyz,
                              double* d2psi_dzz);

__global__ void dot_product_stress(double* d2psi_dxx,
                                   double* d2psi_dxy,
                                   double* d2psi_dxz,
                                   double* d2psi_dyy,
                                   double* d2psi_dyz,
                                   double* d2psi_dzz,
                                   double* psir_ylm_dm,
                                   double* stress_dot,
                                   int elements_num);

__global__ void dot_product_force(double* __restrict__ dpsi_dx,
                                  double* __restrict__ dpsi_dy,
                                  double* __restrict__ dpsi_dz,
                                  double* __restrict__ psir_ylm_dm,
                                  double* force_dot,
                                  int* iat,
                                  int nwmax);

} // namespace GintKernel
#endif // GINT_VL_CUH
