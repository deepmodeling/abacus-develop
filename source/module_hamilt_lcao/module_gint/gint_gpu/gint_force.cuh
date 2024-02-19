#ifndef GINT_FORCE_CUH
#define GINT_FORCE_CUH

#include <cuda_runtime.h>

__global__ void get_psi_force(double *ylmcoef,
                        double delta_r_g,
                        double bxyz_g,
                        double nwmax_g,
                        double *input_double,
                        int *input_int,
                        int *num_psir,
                        int psi_size_max,
                        int *ucell_atom_nwl,
                        bool *atom_iw2_new,
                        int *atom_iw2_ylm,
                        int *atom_iw2_l,
                        int *atom_nw,
                        int nr_max,
                        double *psi_u,
                        double *psir_ylm_right,
                        double *dpsir_ylm_left_x,
                        double *dpsir_ylm_left_y,
                        double *dpsir_ylm_left_z,
                        double *ddpsir_ylm_left_xx,
                        double *ddpsir_ylm_left_xy,
                        double *ddpsir_ylm_left_xz,
                        double *ddpsir_ylm_left_yy,
                        double *ddpsir_ylm_left_yz,
                        double *ddpsir_ylm_left_zz);

__global__ void psir_dot_stress(int * n,
                        double **x_array_g,
                        int incx,
                        double **y_array_g,
                        int incy,
                        double **results_g,
                        int batchcount);

__global__ void dot_product_stress(double *ddpsir_ylm_left_xx,
                                    double *ddpsir_ylm_left_xy,
                                    double *ddpsir_ylm_left_xz,
                                    double *ddpsir_ylm_left_yy,
                                    double *ddpsir_ylm_left_yz,
                                    double *ddpsir_ylm_left_zz,
                                    double *psir_ylm_dm,
                                    double *stress_dot,
                                    int elements_num);

__global__ void dot_product_force(double *dpsir_ylm_left_x,
                                  double *dpsir_ylm_left_y,
                                  double *dpsir_ylm_left_z,
                                  double *psir_ylm_dm,
                                  double *force_dot,
                                  int *iat,
                                  int nwmax,
                                  int max_size,
                                  int elements_num);


#endif // GINT_VL_CUH