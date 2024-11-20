#include "fft_base.h"
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>
#ifndef FFT_ROCM_H
#define FFT_ROCM_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_ROCM : public FFT_BASE<FPTYPE>
{
    public:
        FFT_ROCM(){};
        ~FFT_ROCM(){}; 
        
        void setupFFT() override; 

        void clear() override;

        void cleanFFT() override;

        /**
         * @brief Initialize the fft parameters for ROCM
         * @param nx_in  number of grid points in x direction
         * @param ny_in  number of grid points in y direction
         * @param nz_in  number of grid points in z direction
         * 
         */
        void initfft(int nx_in, 
                    int ny_in, 
                    int nz_in) override;
        
        /**
         * @brief Get the real space data
         * @return real space data
         */
        std::complex<FPTYPE>* get_auxr_3d_data() const override;
        
        /**
         * @brief Forward FFT in 3D for ROCM
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         */
        void fft3D_forward(std::complex<FPTYPE>* in, 
                           std::complex<FPTYPE>* out) const override;

        /**
         * @brief Backward FFT in 3D for ROCM
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         */
        void fft3D_backward(std::complex<FPTYPE>* in, 
                            std::complex<FPTYPE>* out) const override;
    private:
        hipfftHandle c_handle = {};
        hipfftHandle z_handle = {};
        mutable std::complex<float>* c_auxr_3d = nullptr;  // fft space
        mutable std::complex<double>* z_auxr_3d = nullptr; // fft space

};
template FFT_RCOM<float>::FFT_RCOM();
template FFT_ROCM<float>::~FFT_ROCM();
template FFT_RCOM<double>::FFT_RCOM();
template FFT_ROCM<double>::~FFT_ROCM();
}// namespace ModulePW
#endif
