#ifndef GINT_VL_H
#define GINT_VL_H
#include <cublas_v2.h>
#include <cuda.h> // for CUDA_VERSION
#include <cuda_runtime.h>

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

cudaError_t checkCuda(cudaError_t result);

namespace GintKernel
{

void gint_gamma_vl_gpu(hamilt::HContainer<double>* hRGint,
                       const int lgd_now,
                       const int max_size,
                       double vfactor,
                       const double* vlocal,
                       const double* ylmcoef_now,
                       const int pwnczp,
                       const int nbxx,
                       const double dr,
                       const double* rcut,
                       const Grid_Technique& gridt,
                       const UnitCell& ucell);

void gtask_vlocal(const Grid_Technique& gridt,
                  const double* rcut,
                  const UnitCell& ucell,
                  std::vector<bool>& gpu_matrix_calc_flag,
                  const int grid_index_ij,
                  const int max_size,
                  const int nczp,
                  const double vfactor,
                  const double* vlocal_global_value,
                  double* input_double,
                  int* input_int,
                  int* num_psir);

void alloc_mult_vlocal(const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        std::vector<bool>& gpu_matrix_calc_flag,
                        const int grid_index_ij,
                        const int max_size,
                        double* psir_ylm_left,
                        double* psir_r,
                        int* atom_pair_A_m,
                        int* atom_pair_B_n,
                        int* atom_pair_lda,
                        int* atom_pair_ldb,
                        int* atom_pair_ldc,
                        double** atom_pair_mat_A,
                        double** atom_pair_mat_B,
                        double** atom_pair_mat_C,
                        int& atom_pair_num,
                        int& max_m,
                        int& max_n);
} // namespace GintKernel

#endif