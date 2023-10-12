#ifndef GINT_VL_H
#define GINT_VL_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"

void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint, int lgd_now,
                       int nnnmax, const int max_size, double vfactor,
                       const double *vlocal, const double *ylmcoef_now,
                       int pwbx, int pwby, int pwbz, int pwbxyz, int pwncx,
                       int pwncy, int pwnczp, int NLOCAL_now, int nbxx,
                       int *start_ind, const Grid_Technique &GridT);

void gpu_task_generate_vlocal(const Grid_Technique &GridT, int i, int j, int bx,
                              int by, int bz, int bxyz,
                              int atom_pair_size_of_meshcell, int psi_size_max,
                              int max_size, const int ncx, const int ncy,
                              const int nczp, double *dr, int *it,
                              int *psir_ylm_start, int *ib_index, int *num_psir,
                              int *start_ind, int *vindex,
                              int *atom_pair_input_info, int *num_atom_pair);

#endif