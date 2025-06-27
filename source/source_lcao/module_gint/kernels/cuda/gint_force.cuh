#ifndef GINT_FORCE_CUH
#define GINT_FORCE_CUH

#include <cuda_runtime.h>
#include <cstdint>
namespace GintKernel
{

__global__ void get_psi_force(double* ylmcoef,
                              double delta_r,
                              int bxyz,
                              const int nwmax,
                              const int max_atom,
                              const int* const ucell_atom_nwl,
                              const bool* const atom_iw2_new,
                              const int* const atom_iw2_ylm,
                              const int* const atom_iw2_l,
                              const int* const atom_nw,
                              const double* const rcut,
                              const int nr_max,
                              const double* const psi_u,
                              const double* const mcell_pos,
                              const double* const dr_part,
                              const double* const vldr3,
                              const uint8_t* const atoms_type,
                              const int* const atoms_num_info,
                              double* psi,
                              double* dpsi,
                              double* d2psi);


__global__ void dot_product_stress(const double* d2psi,
                                   const double* psi_dm,
                                   const int size,
                                   double* stress);

__global__ void dot_product_force(const int bxyz,
                                  const int nwmax,
                                  const int *atoms_num_info,
                                  const int *iat_on_nbz,
                                  const double* dpsi,
                                  const double* psi_dm,
                                  double* force);

} // namespace GintKernel
#endif // GINT_VL_CUH
