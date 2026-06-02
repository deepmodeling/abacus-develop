#ifndef FFT_BASE_H
#define FFT_BASE_H

#include <complex>
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
     *
     * These virtuals carry a trivial default body so that the vtable always
     * has a valid (non-null) entry for every backend. The previous
     * `__attribute__((weak))` + no-definition pattern relied on the ELF
     * linker resolving unbound weak symbols to null; on Windows/PE (MinGW)
     * that produces null vtable slots and crashes when a non-overridden
     * slot is dispatched. Derived backends still override what they use.
     */
    virtual void initfft(int nx_in,
                         int ny_in,
                         int nz_in,
                         int lixy_in,
                         int rixy_in,
                         int ns_in,
                         int nplane_in,
                         int nproc_in,
                         bool gamma_only_in,
                         bool xprime_in = true) {};

    virtual void initfft(int nx_in, int ny_in, int nz_in) {};

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
    virtual FPTYPE* get_rspace_data() const { return nullptr; }

    virtual std::complex<FPTYPE>* get_auxr_data() const { return nullptr; }

    virtual std::complex<FPTYPE>* get_auxg_data() const { return nullptr; }

    /**
     * @brief Get the auxiliary real space data in 3D
     *
     * The function is used to get the auxiliary real space data in 3D.
     * While the FFT_BASE is an abstract class,the function will be override,
     * The attribute weak is used to avoid define the function.
     */
    virtual std::complex<FPTYPE>* get_auxr_3d_data() const { return nullptr; }

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
    virtual void fftxyfor(std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out) const {};

    virtual void fftxybac(std::complex<FPTYPE>* in,
                          std::complex<FPTYPE>* out) const {};

    /**
     * @brief Forward FFT in z direction
     * @param in  input data
     * @param out  output data
     *
     * This function performs the forward FFT in the z direction.
     * It involves only one axis, z. The FFT is applied only once.
     * Notably, the Z axis operates in many FFT with nz*ns.
     */
    virtual void fftzfor(std::complex<FPTYPE>* in,
                         std::complex<FPTYPE>* out) const {};

    virtual void fftzbac(std::complex<FPTYPE>* in,
                         std::complex<FPTYPE>* out) const {};

    /**
     * @brief Forward FFT in x-y direction with real to complex
     * @param in  input data, real type
     * @param out  output data, complex type
     *
     * This function performs the forward FFT in the x-y direction
     * with real to complex.There is no difference between fftxyfor.
     */
    virtual void fftxyr2c(FPTYPE* in,
                          std::complex<FPTYPE>* out) const {};

    virtual void fftxyc2r(std::complex<FPTYPE>* in,
                          FPTYPE* out) const {};

    /**
     * @brief Forward FFT in 3D
     * @param in  input data
     * @param out  output data
     *
     * This function performs the forward FFT for gpu-like fft.
     * It involves three axes, x, y, and z. The FFT is applied multiple times
     * for fft3D_forward.
     */
    virtual void fft3D_forward(std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out) const {};

    virtual void fft3D_backward(std::complex<FPTYPE>* in,
                                std::complex<FPTYPE>* out) const {};

  protected:
    int nx = 0;
    int ny = 0;
    int nz = 0;
};

template FFT_BASE<float>::FFT_BASE();
template FFT_BASE<double>::FFT_BASE();
template FFT_BASE<float>::~FFT_BASE();
template FFT_BASE<double>::~FFT_BASE();
} // namespace ModuleBase
#endif // FFT_BASE_H
