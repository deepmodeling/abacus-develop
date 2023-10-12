#ifndef GINT_VL_H
#define GINT_VL_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_hamilt_lcao/module_gint/gint.h"

void gint_gamma_vl_gpu(hamilt::HContainer<double> *hRGint, int lgd_now,
                       const int max_size, double vfactor,
                       const double *vlocal, const double *ylmcoef_now,
                       int pwbx, int pwby, int pwbz, int pwbxyz, int pwncx,
                       int pwncy, int pwnczp, int NLOCAL_now, int nbxx,
                       int *start_ind, const Grid_Technique &GridT);

void gpu_task_generate_vlocal(const Grid_Technique &GridT, const int i,
                              const int j, const int bx, const int by,
                              const int bz, const int bxyz,
                              const int atom_pair_size_of_meshcell,
                              const int psi_size_max, const int max_size,
                              const int ncx, const int ncy, const int nczp,
                              const double vfactor,
                              const int *start_ind,
                              const double *vlocal_global_value,
                              double *psi_input_double, int *psi_input_int, int *num_psir,
                              int *atom_pair_input_info, int *num_atom_pair);

#endif