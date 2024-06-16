#ifndef GINT_VL_GPU_H
#define GINT_VL_GPU_H
#include <cublas_v2.h>
#include <cuda.h> // for CUDA_VERSION
#include <cuda_runtime.h>

#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

namespace GintKernel
{

void gint_gamma_vl_gpu(hamilt::HContainer<double>* hRGint,
                       const double* vlocal,
                       const double* ylmcoef_now,
                       const double dr,
                       const double* rcut,
                       const Grid_Technique& gridt,
                       const UnitCell& ucell);

void gtask_vlocal(const Grid_Technique& gridt,
                  const double* rcut,
                  const UnitCell& ucell,
                  std::vector<bool>& gpu_matrix_calc_flag,
                  const int grid_index_ij,
                  const int max_atom,
                  const int nczp,
                  const double vfactor,
                  const double* vlocal_global_value,
                  double* psi_input_double,
                  int* psi_input_int,
                  int* atom_num_per_bcell);

void alloc_mult_vlocal(const Grid_Technique& gridt,
                        const UnitCell& ucell,
                        std::vector<bool>& gpu_matrix_calc_flag,
                        const int grid_index_ij,
                        const int max_atom,
                        double* psi,
                        double* psi_vldr3,
                        std::vector<Cuda_Mem_Wrapper<double>>& grid_vlocal_g,
                        int* mat_m,
                        int* mat_n,
                        int* mat_k,
                        int* mat_lda,
                        int* mat_ldb,
                        int* mat_ldc,
                        double** mat_A,
                        double** mat_B,
                        double** mat_C,
                        int& atom_pair_num,
                        int& max_m,
                        int& max_n);
} // namespace GintKernel

#endif