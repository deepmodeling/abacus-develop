#include "fft_swdfti.h"

#include <cstring>
#include <cstdlib>
#include <mutex>
extern "C" {
#include "swfft.h"   // xMath-SACA swFFT DFTI API (CPE)
}

namespace ModuleBase
{

template <>
void FFT_SWDFTI<double>::setupFFT()
{
    // build all the FFTW plans / buffers exactly as the CPU backend does
    FFT_CPU<double>::setupFFT();

    if (std::getenv("ABACUS_NO_DFTI") != nullptr) { return; }   // A/B: keep FFTW

    // thread-safe one-time CPE spawn (setupFFT may be reached from >1 thread)
    static std::once_flag dfti_athread_once;
    std::call_once(dfti_athread_once, []() { DftiInitAthread(DFTI_SPAWN_QUICK); });

    // batched 1D-z: ns transforms of length nz, contiguous (stride 1, distance nz), in-place
    DFTI_DESCRIPTOR_HANDLE hz = nullptr;
    DftiCreateDescriptor(&hz, DFTI_DOUBLE, DFTI_COMPLEX, 1, (DFTI_LONG)this->nz);
    DftiSetValue(hz, DFTI_NUMBER_OF_TRANSFORMS, (DFTI_LONG)this->ns);
    DftiSetValue(hz, DFTI_INPUT_DISTANCE,  (DFTI_LONG)this->nz);
    DftiSetValue(hz, DFTI_OUTPUT_DISTANCE, (DFTI_LONG)this->nz);
    DftiSetValue(hz, DFTI_PLACEMENT, (DFTI_LONG)DFTI_INPLACE);
    DftiCommitDescriptor(hz);
    this->dftiz = (void*)hz;

    // strided 1D-x: nx-length, (nplane*ny) transforms, stride npy, distance 1
    // (only the xprime / non-gamma k-point layout). y stays on FFTW.
    if (this->xprime && !this->gamma_only)
    {
        const int npy_ = this->nplane * this->ny;
        DFTI_DESCRIPTOR_HANDLE hx = nullptr;
        DftiCreateDescriptor(&hx, DFTI_DOUBLE, DFTI_COMPLEX, 1, (DFTI_LONG)this->nx);
        DftiSetValue(hx, DFTI_NUMBER_OF_TRANSFORMS, (DFTI_LONG)npy_);
        { DFTI_LONG st[2] = {0, (DFTI_LONG)npy_}; DftiSetValue(hx, DFTI_INPUT_STRIDES, st); DftiSetValue(hx, DFTI_OUTPUT_STRIDES, st); }
        DftiSetValue(hx, DFTI_INPUT_DISTANCE,  (DFTI_LONG)1);
        DftiSetValue(hx, DFTI_OUTPUT_DISTANCE, (DFTI_LONG)1);
        DftiSetValue(hx, DFTI_PLACEMENT, (DFTI_LONG)DFTI_INPLACE);
        DftiCommitDescriptor(hx);
        this->dftix = (void*)hx;
    }
}

template <>
void FFT_SWDFTI<double>::cleanFFT()
{
    FFT_CPU<double>::cleanFFT();
    // release the DFTI descriptors before dropping the handles (else they leak)
    if (this->dftiz != nullptr)
    {
        DFTI_DESCRIPTOR_HANDLE hz = (DFTI_DESCRIPTOR_HANDLE)this->dftiz;
        DftiFreeDescriptor(&hz);
        this->dftiz = nullptr;
    }
    if (this->dftix != nullptr)
    {
        DFTI_DESCRIPTOR_HANDLE hx = (DFTI_DESCRIPTOR_HANDLE)this->dftix;
        DftiFreeDescriptor(&hx);
        this->dftix = nullptr;
    }
}

template <>
void FFT_SWDFTI<double>::fftzfor(std::complex<double>* in, std::complex<double>* out) const
{
    if (this->dftiz == nullptr) { FFT_CPU<double>::fftzfor(in, out); return; }
    if (in != out) std::memcpy(out, in, sizeof(std::complex<double>) * (size_t)this->nz * (size_t)this->ns);
    DftiComputeForward((DFTI_DESCRIPTOR_HANDLE)this->dftiz, (void*)out);
}

template <>
void FFT_SWDFTI<double>::fftzbac(std::complex<double>* in, std::complex<double>* out) const
{
    if (this->dftiz == nullptr) { FFT_CPU<double>::fftzbac(in, out); return; }
    if (in != out) std::memcpy(out, in, sizeof(std::complex<double>) * (size_t)this->nz * (size_t)this->ns);
    DftiComputeBackward((DFTI_DESCRIPTOR_HANDLE)this->dftiz, (void*)out);
}

template <>
void FFT_SWDFTI<double>::fftxyfor(std::complex<double>* in, std::complex<double>* out) const
{
    const int npy = this->nplane * this->ny;
    if (this->xprime && this->dftix != nullptr)
    {
        if (in != out) std::memcpy(out, in, sizeof(std::complex<double>) * (size_t)this->nx * (size_t)npy);
        DftiComputeForward((DFTI_DESCRIPTOR_HANDLE)this->dftix, (void*)out);          // x via CPE
        for (int i = 0; i < this->lixy + 1; ++i)                                      // y via FFTW
            fftw_execute_dft(this->planyfor, (fftw_complex*)&out[i * npy], (fftw_complex*)&out[i * npy]);
        for (int i = this->rixy; i < this->nx; ++i)
            fftw_execute_dft(this->planyfor, (fftw_complex*)&out[i * npy], (fftw_complex*)&out[i * npy]);
        return;
    }
    FFT_CPU<double>::fftxyfor(in, out);   // non-xprime / disabled -> FFTW
}

template <>
void FFT_SWDFTI<double>::fftxybac(std::complex<double>* in, std::complex<double>* out) const
{
    const int npy = this->nplane * this->ny;
    if (this->xprime && this->dftix != nullptr)
    {
        if (in != out) std::memcpy(out, in, sizeof(std::complex<double>) * (size_t)this->nx * (size_t)npy);
        for (int i = 0; i < this->lixy + 1; ++i)                                      // y via FFTW
            fftw_execute_dft(this->planybac, (fftw_complex*)&out[i * npy], (fftw_complex*)&out[i * npy]);
        for (int i = this->rixy; i < this->nx; ++i)
            fftw_execute_dft(this->planybac, (fftw_complex*)&out[i * npy], (fftw_complex*)&out[i * npy]);
        DftiComputeBackward((DFTI_DESCRIPTOR_HANDLE)this->dftix, (void*)out);          // x via CPE
        return;
    }
    FFT_CPU<double>::fftxybac(in, out);   // non-xprime / disabled -> FFTW
}

template class FFT_SWDFTI<double>;
} // namespace ModuleBase
