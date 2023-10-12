#ifndef GINT_VL_H
#define GINT_VL_H
#include "module_hamilt_lcao/module_gint/grid_technique.h"

void gint_gamma_vl_gpu(hamilt::HContainer<double>* hRGint,
                       int lgd_now,
                       int nnnmax,
                       const int max_size,
                       double vfactor,
                       const double* vlocal,
                       const double* ylmcoef_now,
                       int pwbx,
                       int pwby,
                       int pwbz,
                       int pwbxyz,
                       int pwncx,
                       int pwncy,
                       int pwnczp,
                       int NLOCAL_now,                       
                       int nbxx,
                       int* start_ind,
                       const Grid_Technique & GridT);

#endif