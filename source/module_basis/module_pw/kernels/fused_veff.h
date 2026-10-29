#ifndef MODULEPW_FUSED_VEFF_H
#define MODULEPW_FUSED_VEFF_H

#include <complex>

namespace ModulePW
{

struct FusedFFT
{
    void (*potential)(int, int, int, const double*, std::complex<double>*, std::complex<double>*, int);
    void (*density)(int, int, int, double, double*, const std::complex<double>*, const std::complex<double>*);
    void (*sticks)(std::complex<double>*, std::complex<double>*, const int*, int, int, int, int, bool, int, int);
};

// Returns generated kernels for a supported length, or nullptr for cuFFT.
const FusedFFT* fused_fft(int length);

// The fused transform is selected only on the measured Volta architecture.
bool fused_veff_supported();

// The two spin components are separated by stride coefficients. Grid storage
// concatenates the retained left and right x slabs, with z varying fastest.
void pack_spinor_grid(int npw, int stride, int plane, int left, int right, const int* index,
                     const std::complex<double>* input, std::complex<double>* up, std::complex<double>* down);

void add_spinor_grid(int npw, int stride, int plane, int left, int right, int nxyz, const int* index,
                    const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output);

// Stick indices address contiguous z lines. The reverse map contains the
// stick index for each x-y pair, or -1 for an absent Fourier line.
void pack_spinor_sticks(int npw, int stride, const int* index, const std::complex<double>* input,
                        std::complex<double>* up, std::complex<double>* down, int bands = 1, int stick_size = 0);
void add_spinor_sticks(int npw, int stride, int nxyz, const int* index,
                       const std::complex<double>* up, const std::complex<double>* down, std::complex<double>* output, int bands = 1, int stick_size = 0);
void transform_stick_slabs(std::complex<double>* sticks, std::complex<double>* slabs,
                           const int* reverse, int nx, int ny, int nz, int left, int right, bool inverse, int bands = 1, int stick_size = 0);

// Applies the inverse x transform, the Hermitian 2x2 spin potential, and the
// forward x transform in place. Normalization is performed when the final
// plane-wave coefficients are gathered.
void fused_veff_x(int nx, int plane, int left, int right, const double* potential,
                  std::complex<double>* up, std::complex<double>* down, int bands = 1);

// Adds the weighted charge and three magnetization components from the
// inverse x transforms. Density channels are contiguous full real grids.
void fused_density_x(int nx, int plane, int left, int right, double weight, double* density,
                     const std::complex<double>* up, const std::complex<double>* down);

} // namespace ModulePW

#endif
