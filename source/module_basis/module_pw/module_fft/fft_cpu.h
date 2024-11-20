#include "fft_base.h"
#include "fftw3.h"

// #ifdef __ENABLE_FLOAT_FFTW

// #endif
// #endif
#ifndef FFT_CPU_H
#define FFT_CPU_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_CPU : public FFT_BASE<FPTYPE>
{
    public:
    FFT_CPU(){};
    FFT_CPU(const int fft_mode_in):fft_mode(fft_mode_in){};
    ~FFT_CPU(){}; 

    /**
     * @brief Initialize the fft parameters.
     * @param nx_in  number of grid points in x direction.
     * @param ny_in  number of grid points in y direction.
     * @param nz_in  number of grid points in z direction.
     * @param lixy_in  the position of the left boundary
     * in the x-y plane.
     * @param rixy_in  the position of the right boundary
     * in the x-y plane.
     * @param ns_in  number of stick whcih is used in the
     * Z direction.
     * @param nplane_in  number of x-y planes.
     * @param nproc_in  number of processors.
     * @param gamma_only_in  whether only gamma point is used.
     * @param xprime_in  whether xprime is used.
     */
    __attribute__((weak)) 
    void initfft(int nx_in, 
                 int ny_in, 
                 int nz_in, 
                 int lixy_in, 
                 int rixy_in, 
                 int ns_in, 
                 int nplane_in, 
                 int nproc_in, 
                 bool gamma_only_in, 
                 bool xprime_in = true) override;
	__attribute__((weak)) 
    void setupFFT() override; 

	// void initplan(const unsigned int& flag = 0);
    __attribute__((weak)) 
    void cleanFFT() override;

    __attribute__((weak)) 
    void clear() override;

    /**
     * @brief Get the real space data the CPU FFT.
     * @return FPTYPE*  the real space data.
     * 
     * the function will return the real space data,
     * which is used in the CPU fft.Use the weak attribute
     * to avoid defining float while without flag ENABLE_FLOAT_FFTW.
     */
    __attribute__((weak)) 
    FPTYPE* get_rspace_data() const override;

    __attribute__((weak)) 
    std::complex<FPTYPE>* get_auxr_data() const override;

    __attribute__((weak)) 
    std::complex<FPTYPE>* get_auxg_data() const override;

    /**
     * @brief Forward FFT in x-y direction
     * @param in  input data
     * @param out  output data
     * 
     * The function details can be found in FFT_BASE,
     * and the function interfaces can be found in FFT_BUNDLE.
     */
    __attribute__((weak)) 
    void fftxyfor(std::complex<FPTYPE>* in, 
                  std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) 
    void fftxybac(std::complex<FPTYPE>* in, 
                  std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) 
    void fftzfor(std::complex<FPTYPE>* in, 
                  std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) 
    void fftzbac(std::complex<FPTYPE>* in, 
                 std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) 
    void fftxyr2c(FPTYPE* in, 
                  std::complex<FPTYPE>* out) const override;

    __attribute__((weak)) 
    void fftxyc2r(std::complex<FPTYPE>* in, 
                  FPTYPE* out) const override;
    private:
        void clearfft(fftw_plan& plan);
        void clearfft(fftwf_plan& plan);

        fftw_plan planzfor  = NULL;
        fftw_plan planzbac  = NULL;
        fftw_plan planxfor1 = NULL;
        fftw_plan planxbac1 = NULL;
        fftw_plan planxfor2 = NULL;
        fftw_plan planxbac2 = NULL;
        fftw_plan planyfor  = NULL;
        fftw_plan planybac  = NULL;
        fftw_plan planxr2c  = NULL;
        fftw_plan planxc2r  = NULL;
        fftw_plan planyr2c  = NULL;
        fftw_plan planyc2r  = NULL;

        fftwf_plan planfzfor = NULL;
        fftwf_plan planfzbac = NULL;
        fftwf_plan planfxfor1= NULL;
        fftwf_plan planfxbac1= NULL;
        fftwf_plan planfxfor2= NULL;
        fftwf_plan planfxbac2= NULL;
        fftwf_plan planfyfor = NULL;
        fftwf_plan planfybac = NULL;
        fftwf_plan planfxr2c = NULL;
        fftwf_plan planfxc2r = NULL;
        fftwf_plan planfyr2c = NULL;
        fftwf_plan planfyc2r = NULL;
        
        std::complex<float>*c_auxg = nullptr;
        std::complex<float>*c_auxr = nullptr;  // fft space
        std::complex<double>*z_auxg = nullptr;
        std::complex<double>*z_auxr = nullptr; // fft space

        float* s_rspace = nullptr;  // real number space for r, [nplane * nx *ny]
        double* d_rspace = nullptr; // real number space for r, [nplane * nx *ny]
        int fftnx=0;
        int fftny=0;
        int fftnxy=0;
        int nxy=0;
        int nplane=0;
        int ns=0; //number of sticks
        int nproc=1; // number of proc.
        int maxgrids = 0;  
        bool gamma_only = false;

        /**
         * @brief lixy: the left edge of the pw ball in the y direction
         */
        int lixy=0;

        /**
         * @brief rixy: the right edge of the pw ball in the x or y direction
         */
        int rixy=0;
        /**
         * @brief xprime: whether xprime is used,when do recip2real, x-fft will 
         * be done last and when doing real2recip, x-fft will be done first; 
         * false: y-fft For gamma_only, true: we use half x; false: we use half y
         */
        bool xprime = true;
        /**
         * @brief fft_mode: fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive
         */
        int fft_mode = 0; 
};
}
#endif // FFT_CPU_H