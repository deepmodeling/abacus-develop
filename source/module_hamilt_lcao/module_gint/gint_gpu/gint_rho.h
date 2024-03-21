#ifndef GINT_RHO_H
#define GINT_RHO_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>

cudaError_t checkCuda(cudaError_t result);
namespace lcaoCudaKernel{

void gint_gamma_rho_gpu(hamilt::HContainer<double> *dm,
                        double *rho,
                        const int nczp,
                        const double *ylmcoef_now,
                        const LCAO_Orbitals &ORB,
                        const Grid_Technique &gridt,
                        const UnitCell &ucell);

void gpu_task_generator_rho(const Grid_Technique &gridt, 
                            const int i, const int j,
                            const int psi_size_max, const int max_size,
                            const int nczp,
                            const UnitCell &ucell,
                            const LCAO_Orbitals &ORB,
                            double *psi_input_double, int *psi_input_int,
                            int *num_psir,
                            const int lgd,
                            double *psir_ylm_g,
                            double *psir_zeros_g,
                            double *dm_matrix_g,
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

void rho_cal_task( const Grid_Technique &gridt,
                      const int i, const int j,
                      const int max_size,
                      const int lgd,
                      const bool *gpu_mat_cal_flag,
                      const int *start_idx_psi,
                      const double *psir_ylm,
                      const double *psir_zeros,
                      const double *dm_matrix,
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
                      int &atom_pair_num
                      );

} // namespace lcaoCudaKernel
#endif