#ifndef GINT_VL_H
#define GINT_VL_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/gint_gpu/vbatch_matrix_multiple/cuda_tools.cuh"

void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint, int lgd_now,
                       const int max_size, double vfactor,
                       const double *vlocal, const double *ylmcoef_now,
                       int pwnczp, int NLOCAL_now, int nbxx, const Grid_Technique &GridT);




void gpu_task_generate_vlocal(const Grid_Technique &GridT, 
                              const int i, const int j,
                              const int atom_pair_size_of_meshcell_v2,
                              const int psi_size_max, const int max_size,
                              const int nczp,
                              const double vfactor,
                              const double *vlocal_global_value,
                              double *psir_ylm_left,
                              double *psir_ylm_right,
                              double *psi_input_double, int *psi_input_int,
                              int *num_psir, int *atom_pair_input_info,
                              int *num_atom_pair, double* GridVlocal_v2_g[],     
                              double ** atom_pair_left_v2,
                              double ** atom_pair_right_v2,
                              double ** atom_pair_output_v2);

#endif