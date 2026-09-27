#ifndef ABACUS_GREENX_MINIMAX_ADAPTER_H
#define ABACUS_GREENX_MINIMAX_ADAPTER_H

#ifdef __GREENX_MINIMAX
extern "C"
{
void gx_minimax_grid_frequency_wrp(int* num_points,
                                   double* e_min,
                                   double* e_max,
                                   double* omega_points,
                                   double* omega_weights,
                                   int* ierr);
}
#endif

#endif
