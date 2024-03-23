#ifndef GINT_RHO_H
#define GINT_RHO_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>

cudaError_t checkCuda(cudaError_t result);
namespace lcaoCudaKernel{

/**
 * calculate the rho by GPU
 *
 * @param dm density matrix.
 * @param nczp number of meshcells along the z-axis on this processor.
 * @param ylmcoef_now coefficients for the spherical harmonics expansion.
 * @param gridt Grid_Technique object containing grid information.
 * @param ORB LCAO_Orbitals.
 * @param ucell UnitCell.
 * @param rho array to store rho.
 */
void gint_gamma_rho_gpu(const hamilt::HContainer<double> *dm,
                        const int nczp,
                        const double *ylmcoef_now,
                        const Grid_Technique &gridt,
                        const LCAO_Orbitals &ORB,
                        const UnitCell &ucell,
                        double *rho);

/**
 * generate GPU tasks for computing the rho.
 * the computation task can be divided into psir calculation, matrix multiplication and vector dot product.
 * the matrix multiplication is mat_dm * mat_psir, and the vector dot product is psir * psir_dm.
 * This function will be split into three separate functions,
 * which are calculating psir, matrix multiplication, and vector dot product.
 *
 * @param gridt Grid_Technique object containing grid information.
 * @param i X index of the bigcell.
 * @param j Y index of the bigcell.
 * @param max_size maximum number of atoms on a meshcell.
 * @param nczp number of meshcells along the z-axis on this processor.
 * @param ucell UnitCell object containing unit cell information.
 * @param ORB LCAO_Orbitals object containing LCAO orbital information.
 * @param psi_input_double array storing `double` type data used for calculating psir.
 * @param psi_input_int array storing `int` type data used for calculating psir.
 * @param num_psir array records the number of atoms on each bigcell.
 * @param lgd lgd.
 * @param psir_ylm_g array used to store psir.
 * @param psir_dm_g array used to store psir_dm.
 * @param dm_matrix_g array used to store  dm_matrix.
 * @param mat_alpha array containing alpha values for matrix multiplication.
 * @param mat_m array containing m values for matrix multiplication.
 * @param mat_n array containing n values for matrix multiplication.
 * @param mat_k array containing k values for matrix multiplication.
 * @param mat_lda array containing lda values for matrix multiplication.
 * @param mat_ldb array containing ldb values for matrix multiplication.
 * @param mat_ldc array containing ldc values for matrix multiplication.
 * @param mat_A matrix A for matrix multiplication.
 * @param mat_B matrix B for matrix multiplication.
 * @param mat_C matrix C for matrix multiplication.
 * @param max_m maximum value of m.
 * @param max_n maximum value of n.
 * @param atom_pair_num total number of atom pairs, 
 *                      which is also the number of mat mul operations.
 * @param rho_g array used to store rho.
 * @param vec_l psir_ylm for vec dot product.
 * @param vec_r psir_dm for vec dot product.
 * @param dot_product array storing the positions where each dot product should be stored.
 * @param vec_len Array storing the vector lengths for each dot product.
 * @param dot_count total count of dot products.
 */
void gpu_task_generator_rho(const Grid_Technique &gridt, 
                            const int i, const int j,
                            const int max_size,
                            const int nczp,
                            const UnitCell &ucell,
                            const LCAO_Orbitals &ORB,
                            double *psi_input_double, int *psi_input_int,
                            int *num_psir,
                            const int lgd,
                            double * const psir_ylm_g,
                            double * const psir_dm_g,
                            double * const dm_matrix_g,
                            double *mat_alpha,
                            int *mat_m,
                            int *mat_n,
                            int *mat_k,
                            int *mat_lda,
                            int *mat_ldb,
                            int *mat_ldc,
                            double **mat_A,
                            double **mat_B,
                            double **mat_C,
                            int &max_m,
                            int &max_n,
                            int &atom_pair_num,
                            double *rho_g,
                            double **vec_l,
                            double **vec_r,
                            double **dot_product,
                            int *vec_len,
                            int &dot_count 
                            );

} // namespace lcaoCudaKernel
#endif