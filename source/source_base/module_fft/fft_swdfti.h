#ifndef FFT_SWDFTI_H
#define FFT_SWDFTI_H
// CPE-accelerated CPU FFT backend: subclasses FFT_CPU and overrides only the
// local 1D sticks FFTs (batched z, strided x) with the Sunway swFFT xMath DFTI
// API (offloaded to the 64 CPEs via DftiInitAthread). Everything else (plan
// setup, 2D-xy y-direction, r2c/c2r, box 3D) is inherited from FFT_CPU/FFTW.
// Compiled only on Sunway (USE_SWDFTI) and selected by the FFT factory in
// FFT_Bundle for device "cpu" -- so FFT_CPU itself stays free of any DFTI macro.
#include "fft_cpu.h"

namespace ModuleBase
{
template <typename FPTYPE>
class FFT_SWDFTI : public FFT_CPU<FPTYPE>
{
  public:
    FFT_SWDFTI() {};
    FFT_SWDFTI(const int fft_mode_in) : FFT_CPU<FPTYPE>(fft_mode_in) {};
    ~FFT_SWDFTI() {};

    void setupFFT() override;
    void cleanFFT() override;

    void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;
    void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;
    void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;
    void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const override;

  private:
    // swFFT DFTI descriptors: z (batched ns x nz contiguous) and x (strided).
    // y stays on FFTW (CPE loses on the small per-slice y-batch). null => FFTW.
    void* dftiz = nullptr;
    void* dftix = nullptr;
};
} // namespace ModuleBase
#endif
