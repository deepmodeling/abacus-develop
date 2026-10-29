#ifndef FFT_CUDA_H
#define FFT_CUDA_H

#include "fft_base.h"
#include "cufft.h"
#include "cuda_runtime.h"
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

        // Transform only the retained x slabs when the plane-wave support is sparse.
        void fft3D_forward_pw(std::complex<FPTYPE>* in,
                              std::complex<FPTYPE>* out) const;
        void fft3D_backward_pw(std::complex<FPTYPE>* in,
                               std::complex<FPTYPE>* out) const;

#if defined(__CUFFTDX)
        int fused_veff_grid_size() const;
        int fused_veff_stick_size() const;
        std::complex<double>* get_spinor_stick_data(int bands) const;
        void transform_spinor_sticks(std::complex<double>* up, std::complex<double>* down,
                                     const int* reverse, bool inverse, int bands = 1) const;
        void apply_spinor_potential(const double* potential, std::complex<double>* up,
                                    std::complex<double>* down, const int* reverse = nullptr, int bands = 1) const;
#endif
        int lixy() const { return this->lixy_index; }
        int rixy() const { return this->rixy_index; }
    private:
        mutable cufftHandle c_handle = {};
        mutable cufftHandle z_handle = {};
        mutable cufftHandle z_x_handle = {};
        mutable cufftHandle z_yz_left_handle = {};
        mutable cufftHandle z_yz_right_handle = {};
        mutable cufftHandle z_stick_handle[3] = {};
       
        mutable std::complex<float>* c_auxr_3d = nullptr;  // fft space
        mutable std::complex<double>* z_auxr_3d = nullptr; // fft space
        mutable int z_auxr_3d_capacity = 0;

        int lixy_index = 0;
        int rixy_index = 0;
        int nproc = 1;
        int ns_index = 0;
        bool gamma_only = false;
        bool xprime = true;

        bool use_pruned_gpu_fft() const;
        void init_gpu_pruned_plan() const;

};

} // namespace ModulePW
#endif