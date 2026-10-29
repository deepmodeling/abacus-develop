#ifndef MODULEPW_FUSED_VEFF_H
#define MODULEPW_FUSED_VEFF_H

#include <complex>

namespace ModulePW
{

// The fused transform is selected only on the measured Volta architecture.
bool fused_veff_supported();

// The two spin components are separated by stride coefficients. Grid storage
// concatenates the retained left and right x slabs, with z varying fastest.
void pack_spinor_grid(int npw, int stride, int plane, int left, int right, const int* index,
                     const std::complex<double>* input, std::complex<double>* up, std::complex<double>* down);

void add_spinor_grid(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                    const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output);

// Applies the inverse x transform, the Hermitian 2x2 spin potential, and the
// forward x transform in place. The x dimension is 200; normalization is
// performed when the final plane-wave coefficients are gathered.
void fused_veff_x(int plane, int left, int right, const double* potential,
                  std::complex<double>* up, std::complex<double>* down);

} // namespace ModulePW

#endif
