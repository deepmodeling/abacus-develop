#include "fft_base.h"
#include "cufft.h"
#include "cuda_runtime.h"

#ifndef FFT_CUDA_H
#define FFT_CUDA_H
namespace ModulePW
{
template <typename FPTYPE>
class FFT_CUDA : public FFT_BASE<FPTYPE>
{
    public:
        FFT_CUDA(){};
        ~FFT_CUDA(){}; 
        
	    void setupFFT() override; 

        void clear() override;

        void cleanFFT() override;

        /** 
        * @brief Initialize the fft parameters
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
         * @brief Forward FFT in 3D
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         * 
         * This function performs the forward FFT in 3D.
         */
        void fft3D_forward(std::complex<FPTYPE>* in, 
                           std::complex<FPTYPE>* out) const override;
        /**
         * @brief Backward FFT in 3D
         * @param in  input data, complex FPTYPE
         * @param out  output data, complex FPTYPE
         * 
         * This function performs the backward FFT in 3D.
         */
        void fft3D_backward(std::complex<FPTYPE>* in, 
                            std::complex<FPTYPE>* out) const override;
    private:
        cufftHandle c_handle = {};
        cufftHandle z_handle = {};
       
        std::complex<float>* c_auxr_3d = nullptr;  // fft space
        std::complex<double>* z_auxr_3d = nullptr; // fft space

};
template FFT_CUDA<float>::FFT_CUDA();
template FFT_CUDA<float>::~FFT_CUDA();
template FFT_CUDA<double>::FFT_CUDA();
template FFT_CUDA<double>::~FFT_CUDA();
} // namespace ModulePW
#endif