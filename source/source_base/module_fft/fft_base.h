#ifndef FFT_BASE_H
#define FFT_BASE_H

#include <complex>

// These FFT virtuals are declared weak so the ELF linker can resolve the
// unused single-precision (FFT_CPU<float>) vtable slots to null when
// ENABLE_FLOAT_FFTW is off. MinGW/PE has no working equivalent: weak template
// members there either collide ("multiple definition") or leave null vtable
// slots that crash on dispatch. On Windows we therefore drop the attribute and
// rely on the build defining the symbols (ENABLE_FLOAT_FFTW=ON supplies the
// real FFT_CPU<float> methods; the float CPU path is unused otherwise).
// Linux/ELF behaviour is unchanged -- ABACUS_FFT_WEAK expands to exactly
// __attribute__((weak)) there.
#if defined(_WIN32)
#define ABACUS_FFT_WEAK
#else
#define ABACUS_FFT_WEAK __attribute__((weak))
#endif

namespace ModuleBase
{
template <typename FPTYPE>
class FFT_BASE
{
  public:
    FFT_BASE() {};
    virtual ~FFT_BASE() {};

    /**
     * @brief Initialize the fft parameters as virtual function.
     *
     * The function is used to initialize the fft parameters.
     */
    virtual ABACUS_FFT_WEAK void initfft(int nx_in,
                                               int ny_in,
                                               int nz_in,
                                               int lixy_in,
                                               int rixy_in,
                                               int ns_in,
                                               int nplane_in,
                                               int nproc_in,
                                               bool gamma_only_in,
                                               bool xprime_in = true);

    virtual ABACUS_FFT_WEAK void initfft(int nx_in, int ny_in, int nz_in);

    /**
     * @brief Setup the fft plan and data as pure virtual function.
     *
     * The function is set as pure virtual function.In order to
     * override the function in the derived class.In the derived
     * class, the function is used to setup the fft plan and data.
     */
    virtual void setupFFT() = 0;

    /**
     * @brief Clean the fft plan as pure virtual function.
     *
     * The function is set as pure virtual function.In order to
     * override the function in the derived class.In the derived
     * class, the function is used to clean the fft plan.
     */
    virtual void cleanFFT() = 0;

    /**
     * @brief Clear the fft data as pure virtual function.
     *
     * The function is set as pure virtual function.In order to
     * override the function in the derived class.In the derived
     * class, the function is used to clear the fft data.
     */
    virtual void clear() = 0;
    /**
     * @brief Allocate and destory the resoure in FFT running time,
     * Now it only used in the DSP mode.
     * 
     * The function is set as pure virtual function.In order to
     * override the function in the derived class.In the derived
     * class, the function is used to allocate and destory the
     * resoure in FFT running time.
     */
    virtual void resource_handler(const int flag) const {};
    /**
     * @brief Get the real space data in cpu-like fft
     *
     * The function is used to get the real space data.While the
     * FFT_BASE is an abstract class,the function will be override,
     * The attribute weak is used to avoid define the function.
     */
    virtual ABACUS_FFT_WEAK FPTYPE* get_rspace_data() const;

    virtual ABACUS_FFT_WEAK std::complex<FPTYPE>* get_auxr_data() const;

    virtual ABACUS_FFT_WEAK std::complex<FPTYPE>* get_auxg_data() const;

    /**
     * @brief Get the auxiliary real space data in 3D
     *
     * The function is used to get the auxiliary real space data in 3D.
     * While the FFT_BASE is an abstract class,the function will be override,
     * The attribute weak is used to avoid define the function.
     */
    virtual ABACUS_FFT_WEAK std::complex<FPTYPE>* get_auxr_3d_data() const;

    // forward fft in x-y direction

    /**
     * @brief Forward FFT in x-y direction
     * @param in  input data
     * @param out  output data
     *
     * This function performs the forward FFT in the x-y direction.
     * It involves two axes, x and y. The FFT is applied multiple times
     * along the left and right boundaries in the primary direction(which is
     * determined by the xprime flag).Notably, the Y axis operates in
     * "many-many-FFT" mode.
     */
    virtual ABACUS_FFT_WEAK void fftxyfor(std::complex<FPTYPE>* in, 
                                                std::complex<FPTYPE>* out) const;

    virtual ABACUS_FFT_WEAK void fftxybac(std::complex<FPTYPE>* in, 
                                                std::complex<FPTYPE>* out) const;

    /**
     * @brief Forward FFT in z direction
     * @param in  input data
     * @param out  output data
     *
     * This function performs the forward FFT in the z direction.
     * It involves only one axis, z. The FFT is applied only once.
     * Notably, the Z axis operates in many FFT with nz*ns.
     */
    virtual ABACUS_FFT_WEAK void fftzfor(std::complex<FPTYPE>* in, 
                                               std::complex<FPTYPE>* out) const;

    virtual ABACUS_FFT_WEAK void fftzbac(std::complex<FPTYPE>* in, 
                                               std::complex<FPTYPE>* out) const;

    /**
     * @brief Forward FFT in x-y direction with real to complex
     * @param in  input data, real type
     * @param out  output data, complex type
     *
     * This function performs the forward FFT in the x-y direction
     * with real to complex.There is no difference between fftxyfor.
     */
    virtual ABACUS_FFT_WEAK void fftxyr2c(FPTYPE* in, 
                                                std::complex<FPTYPE>* out) const;

    virtual ABACUS_FFT_WEAK void fftxyc2r(std::complex<FPTYPE>* in, 
                                                FPTYPE* out) const;

    /**
     * @brief Forward FFT in 3D
     * @param in  input data
     * @param out  output data
     *
     * This function performs the forward FFT for gpu-like fft.
     * It involves three axes, x, y, and z. The FFT is applied multiple times
     * for fft3D_forward.
     */
    virtual ABACUS_FFT_WEAK void fft3D_forward(std::complex<FPTYPE>* in, 
                                                     std::complex<FPTYPE>* out) const;

    virtual ABACUS_FFT_WEAK void fft3D_backward(std::complex<FPTYPE>* in, 
                                                      std::complex<FPTYPE>* out) const;

  protected:
    int nx = 0;
    int ny = 0;
    int nz = 0;
};

#if defined(_WIN32)
// On Linux the non-pure base virtuals above are __attribute__((weak)) and the
// ELF linker resolves their (never-used) vtable slots to null. MinGW/PE has no
// such fallback, so define trivial bodies for them here -- they are never
// executed (FFT_BASE is abstract; every backend overrides what it actually
// uses, and the unoverridden slots, e.g. fft3D_* on the CPU backend, are not
// called). This block is compiled only on Windows; Linux keeps the upstream
// weak declarations unchanged.
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::initfft(int, int, int, int, int, int, int, int, bool, bool) {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::initfft(int, int, int) {}
template <typename FPTYPE>
FPTYPE* FFT_BASE<FPTYPE>::get_rspace_data() const { return nullptr; }
template <typename FPTYPE>
std::complex<FPTYPE>* FFT_BASE<FPTYPE>::get_auxr_data() const { return nullptr; }
template <typename FPTYPE>
std::complex<FPTYPE>* FFT_BASE<FPTYPE>::get_auxg_data() const { return nullptr; }
template <typename FPTYPE>
std::complex<FPTYPE>* FFT_BASE<FPTYPE>::get_auxr_3d_data() const { return nullptr; }
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftxyfor(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftxybac(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftzfor(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftzbac(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftxyr2c(FPTYPE*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fftxyc2r(std::complex<FPTYPE>*, FPTYPE*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fft3D_forward(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
template <typename FPTYPE>
void FFT_BASE<FPTYPE>::fft3D_backward(std::complex<FPTYPE>*, std::complex<FPTYPE>*) const {}
#endif // _WIN32

template FFT_BASE<float>::FFT_BASE();
template FFT_BASE<double>::FFT_BASE();
template FFT_BASE<float>::~FFT_BASE();
template FFT_BASE<double>::~FFT_BASE();
} // namespace ModuleBase
#endif // FFT_BASE_H
