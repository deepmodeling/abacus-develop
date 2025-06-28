#ifndef GINT_VL_GPU_H
#define GINT_VL_GPU_H

#include "gint.h"
#include "grid_technique.h"
#include "kernels/cuda/cuda_tools.cuh"

namespace GintKernel
{

void gint_vl_gpu(hamilt::HContainer<double>* hRGint,
                 const double* vlocal,
                 const double* ylmcoef_now,
                 const double dr,
                 const double* rcut,
                 const Grid_Technique& gridt,
                 const UnitCell& ucell);

void gtask_vlocal(const Grid_Technique& gridt,
                  const UnitCell& ucell,
                  const int grid_index_ij,
                  const int nczp,
                  const double vfactor,
                  const double* vlocal_global_value,
                  int& atoms_per_z,
                  int* atoms_num_info,
                  uint8_t* atoms_type,
                  double* dr_part,
                  double* vldr3);

void alloc_mult_vlocal(const hamilt::HContainer<double>* hRGint,
                       const Grid_Technique& gridt,
                       const UnitCell& ucell,
                       const int grid_index_ij,
                       const int max_atom,
                       double* const psi,
                       double* const psi_vldr3,
                       double* const grid_vlocal_g,
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