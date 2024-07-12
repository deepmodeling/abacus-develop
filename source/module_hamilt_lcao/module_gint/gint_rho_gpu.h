#ifndef GINT_RHO_GPU_H
#define GINT_RHO_GPU_H
#include "gint.h"

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
void gint_gamma_rho_gpu(const hamilt::HContainer<double>* dm,
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

/**
 * Allocate resources and perform matrix multiplication and vector dot products 
 * for computing the rho.
 *
 * @param gridt Grid_Technique object containing grid information.
 * @param ucell UnitCell object containing unit cell information.
 * @param gpu_mat_cal_flag Flags indicating which parts of the calculation will use the GPU.
 * @param grid_index_ij Combined X and Y indices of the bigcell.
 * @param max_size Maximum number of atoms in a meshcell.
 * @param lgd lgd.
 * @param nczp Number of meshcells along the z-axis on this processor.
 * @param psir_ylm_g One-dimensional array storing psir.
 * @param psir_dm_g One-dimensional array storing psir_dm.
 * @param dm_matrix_g One-dimensional array storing mat_dm.
 * @param mat_alpha Alpha values for matrix multiplication.
 * @param mat_m Number of rows in mat_dm.
 * @param mat_n Number of columns in mat_psir.
 * @param mat_k Number of columns in mat_dm, which equals the number of rows in mat_psir.
 * @param mat_lda Leading dimension of mat_dm.
 * @param mat_ldb Leading dimension of mat_psir.
 * @param mat_ldc Leading dimension of mat_psir_dm.
 * @param mat_A Pointers to mat_dm.
 * @param mat_B Pointers to mat_psir.
 * @param mat_C Pointers to mat_psir_dm.
 * @param max_m Maximum value of m.
 * @param max_n Maximum value of n.
 * @param atom_pair_num Total count of atom pairs, which is also the number of matrix multiplication operations.
 * @param rho_g Rho.
 * @param dot_product Pointers to the results of dot products.
 */
void alloc_mult_dot_rho(const Grid_Technique& gridt,
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