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
void gint_gamma_rho_gpu(const hamilt::HContainer<double>* dm,
                        const double* ylmcoef_now,
                        const double dr,
                        const double* rcut,
                        const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        double* rho);

/**
 * Generate GPU tasks for computing the rho.
 * The computation task can be divided into psir calculation, matrix
 * multiplication, and vector dot product. The matrix multiplication is mat_dm *
 * mat_psir, and the vector dot product is psir * psir_dm. This function will be
 * split into three separate functions: calculating psir, matrix multiplication,
 * and vector dot product.
 *
 * @param gridt Grid_Technique object containing grid information.
 * @param grid_index_ij Combined X and Y indices of the bigcell.
 * @param gpu_matrix_calc_flag Flags indicating which parts of the calculation will use the GPU.
 * @param max_atom Maximum number of atoms in a meshcell.
 * @param ucell UnitCell object containing unit cell information.
 * @param rcut Pointer to the cutoff radius array.
 * @param psi_input_double `double` type data used for calculating psir.
 * @param psi_input_int `int` type data used for calculating psir.
 * @param atom_num_per_bcell Number of atoms in each bigcell.
 */
void gtask_rho(const Grid_Technique& gridt,
               const int grid_index_ij,
               std::vector<bool>& gpu_mat_cal_flag,
               const int max_atom,
               const UnitCell& ucell,
               const double* rcut,
               double* psi_input_double,
               int* psi_input_int,
               int* atom_num_per_bcell,
               int* start_idx_per_bcell,
               int& atom_per_z);

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
                        const std::vector<bool>& gpu_mat_cal_flag,
                        const int grid_index_ij,
                        const int max_size,
                        const int lgd,
                        const int nczp,
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