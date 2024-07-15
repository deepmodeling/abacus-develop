#ifndef GINT_VL_CUH
#define GINT_VL_CUH

#include <cuda_runtime.h>
#include <cstdint>
namespace GintKernel
{
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
__global__ void get_psi_and_vldr3(const double* const ylmcoef,
                                  const double delta_r,
                                  const int bxyz,
                                  const double nwmax,
                                  const double max_atom,
                                  const int* const ucell_atom_nwl,
                                  const bool* const atom_iw2_new,
                                  const int* const atom_iw2_ylm,
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
                                  double* psi_vldr3);

} // namespace GintKernel
#endif // GINT_VL_CUH