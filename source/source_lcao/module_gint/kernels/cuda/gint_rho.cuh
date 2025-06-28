#ifndef GINT_RHO_CUH
#define GINT_RHO_CUH

#include <cuda_runtime.h>
#include <cstdint>
namespace GintKernel
{

/**
 * @brief CUDA kernel to calculate psir.
 *
 * This kernel calculates the wave function psi using the provided input
 * parameters.
 */
__global__ void get_psi(const double* const ylmcoef,
                        const double delta_r,
                        const int bxyz,
                        const int nwmax,
                        const int max_atom,
                        const int* const ucell_atom_nwl,
                        const bool* const atom_iw2_new,
                        const int* const atom_iw2_ylm,
                        const int* const atom_nw,
                        const double* const rcut,
                        const int nr_max,
                        const double* const psi_u,
                        const double* const mcell_pos,
                        const double* const dr_part,
                        const uint8_t* const atoms_type,
                        const int* const atoms_num_info,
                        double* psi);

__global__ void psir_dot(const int bxyz,
                         const int nwmax,
                         const int* atoms_num_info,
                         const double* __restrict__ vec_a_g,
                         const double* __restrict__  vec_b_g,
                         double** results_g);

} // namespace GintKernel
#endif // GINT_RHO_CUH