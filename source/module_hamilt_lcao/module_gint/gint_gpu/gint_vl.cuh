#ifndef GINT_VL_CUH
#define GINT_VL_CUH

#include <cuda_runtime.h>
void gint_gamma_vl_upload_const(const int max_size,
                       const double *ylmcoef_now,
                       const int bxyz);

__global__ void get_psi_and_vldr3(double *input_double,
                                  int *input_int,
                                  int *num_psir,
                                  int psi_size_max,
                                  int *ucell_atom_nwl,
                                  bool *atom_iw2_new,
                                  int *atom_iw2_ylm,
                                  int *atom_nw,
                                  int nr_max,
                                  double *psi_u,
                                  double *psir_ylm_left,
                                  double *psir_ylm_right);

__global__ void psi_multiple(double ** atom_pair_left_g_v2,
                             double ** atom_pair_right_g_v2,
                             double ** atom_pair_output,
                             int *atom_pair_input_info_g,
                             int *num_atom_pair_g,
                             int atom_pair_size_of_meshcell_v2,
                             double *GridVlocal,
                             int lgd);

#endif // GINT_VL_CUH