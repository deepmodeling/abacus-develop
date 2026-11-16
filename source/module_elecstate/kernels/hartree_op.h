#ifndef MODULE_ELECSTATE_HARTREE_OP_H
#define MODULE_ELECSTATE_HARTREE_OP_H

#include <complex>

namespace elecstate
{

int hartree_reduction_size_gpu(int npw);

// Density channels are contiguous real-space fields; their sum is transformed.
void hartree_pack_density_gpu(int size, int channels, const double* density, std::complex<double>* grid);

// The reciprocal scale is evaluated by the host to match PW_Basis normalization.
// Each block writes one partial sum of |rho(G)|^2 times the Coulomb kernel.
void hartree_coefficients_gpu(int npw, const int* index, const double* gg, double scale,
                              double coulomb_numerator, double tpiba2, const std::complex<double>* grid,
                              std::complex<double>* potential, double* energy);

} // namespace elecstate

#endif
