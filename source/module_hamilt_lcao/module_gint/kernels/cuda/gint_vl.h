#ifndef GINT_VL_H
#define GINT_VL_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include <cuda.h>    // for CUDA_VERSION
#include <cublas_v2.h>
#include <cuda_runtime.h>

cudaError_t checkCuda(cudaError_t result);

namespace lcaoCudaKernel{

void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint, int lgd_now,
                       const int max_size, double vfactor,
                       const double *vlocal, const double *ylmcoef_now,
                       int pwnczp, int nbxx, const Grid_Technique &gridt,
                       const LCAO_Orbitals &ORB, const UnitCell &ucell);

void gtask_vlocal(const Grid_Technique &gridt, 
                              const LCAO_Orbitals &ORB,
                              const UnitCell &ucell,
                              const int i, const int j,
                              const int max_size,
                              const int nczp,
                              const double vfactor,
                              const double *vlocal_global_value,
                              double *psir_ylm_left,
                              double *psir_ylm_right,
                              double *psi_input_double, int *psi_input_int,
                              int *num_psir,
                              int *atom_pair_left_info,
                              int *atom_pair_right_info,
                              int *atom_pair_lda,
                              int *atom_pair_ldb,
                              int *atom_pair_ldc,
                              double ** atom_pair_left_v2,
                              double ** atom_pair_right_v2,
                              double ** atom_pair_output_v2,
                              int & atom_pair_num,
                              int & max_m,
                              int & max_n);

} // namespace lcaoCudaKernel

#endif