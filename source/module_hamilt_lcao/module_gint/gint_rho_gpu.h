#ifndef GINT_RHO_H
#define GINT_RHO_H
#include <cublas_v2.h>
#include <cuda.h> // for CUDA_VERSION
#include <cuda_runtime.h>

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

namespace GintKernel
{

/**
 * calculate the rho by GPU
 *
 * @param dm density matrix.
 * @param ylmcoef_now coefficients for the spherical harmonics expansion.
 * @param dr The grid spacing.
 * @param rcut Pointer to the cutoff radius array.
 * @param gridt Grid_Technique object containing grid information.
 * @param ucell UnitCell.
 * @param rho rho.
 */
void gint_rho_gpu(const hamilt::HContainer<double>* dm,
                        const double* ylmcoef_now,
                        const double dr,
                        const double* rcut,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        double* rho);

void gtask_rho(const Grid_Technique& gridt,
               const int grid_index_ij,
               const UnitCell& ucell,
               double* dr_part,
               uint8_t* atoms_type,
               int* atoms_num_info,
               int& atoms_per_z);

void alloc_mult_dot_rho(const hamilt::HContainer<double>* dm,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        const int grid_index_ij,
                        const int max_atom,
                        const int lgd,
                        const int nczp,
                        const int* atoms_num_info,
                        double* const psir_ylm_g,
                        double* const psir_dm_g,
                        double* const dm_matrix_g,
                        double* mat_alpha,
                        int* mat_m,
                        int* mat_n,
                        int* mat_k,
                        int* mat_lda,
                        int* mat_ldb,
                        int* mat_ldc,
                        double** mat_A,
                        double** mat_B,
                        double** mat_C,
                        int& max_m,
                        int& max_n,
                        int& atom_pair_num,
                        double* rho_g,
                        double** dot_product);

} // namespace GintKernel
#endif