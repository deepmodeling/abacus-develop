#ifndef FFT_CUDA_H
#define FFT_CUDA_H

#include "fft_base.h"
#include "cufft.h"
#include "cuda_runtime.h"
namespace ModuleBase
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

        // Batch FFT methods
        /**
         * @brief Setup batch FFT plans and allocate batch buffers
         * @param batch_size_in Number of FFTs per batch (1-128)
         *
         * Must be called after initfft(). Creates cuFFTPlanMany with
         * the specified batch size and allocates device memory.
         */
        void setupBatchFFT(int batch_size_in = 8);

        /**
         * @brief Clean up batch FFT plans
         */
        void cleanBatchFFT();

        /**
         * @brief Forward batch FFT in 3D
         * @param in_batch  input data batch, size = batch_size * nx * ny * nz
         * @param out_batch output data batch, size = batch_size * nx * ny * nz
         * @param batch_count actual number of FFTs to process (may be < BATCH_FFT_SIZE)
         *
         * Performs batch_count forward 3D FFTs in a single kernel launch.
         */
        void fft3D_forward_batch(std::complex<FPTYPE>* in_batch,
                                 std::complex<FPTYPE>* out_batch,
                                 int batch_count) const;

        /**
         * @brief Backward batch FFT in 3D
         * @param in_batch  input data batch, size = batch_size * nx * ny * nz
         * @param out_batch output data batch, size = batch_size * nx * ny * nz
         * @param batch_count actual number of FFTs to process (may be < BATCH_FFT_SIZE)
         *
         * Performs batch_count backward 3D FFTs in a single kernel launch.
         */
        void fft3D_backward_batch(std::complex<FPTYPE>* in_batch,
                                  std::complex<FPTYPE>* out_batch,
                                  int batch_count) const;

        /**
         * @brief Check if batch FFT is ready
         * @return true if batch FFT plans and buffers are initialized
         */
        bool is_batch_fft_ready() const;

        /**
         * @brief Get batch FFT size
         * @return maximum batch size (BATCH_FFT_SIZE)
         */
        int get_batch_size() const;

        /**
         * @brief Get batch input buffer
         * @return pointer to device batch input buffer
         */
        std::complex<FPTYPE>* get_batch_input_buffer() const;

        /**
         * @brief Get batch output buffer
         * @return pointer to device batch output buffer
         */
        std::complex<FPTYPE>* get_batch_output_buffer() const;

    private:
        cufftHandle c_handle = {};
        cufftHandle z_handle = {};

        int batch_size = 8;  // Runtime batch size (configurable)

        std::complex<float>* c_auxr_3d = nullptr;  // fft space
        std::complex<double>* z_auxr_3d = nullptr; // fft space

        // Batch FFT handles and buffers
        cufftHandle c_batch_handle = {};
        cufftHandle z_batch_handle = {};

        std::complex<float>* c_auxr_batch_in = nullptr;   // batch input buffer
        std::complex<float>* c_auxr_batch_out = nullptr;  // batch output buffer
        std::complex<double>* z_auxr_batch_in = nullptr;  // batch input buffer
        std::complex<double>* z_auxr_batch_out = nullptr; // batch output buffer

};

} // namespace ModuleBase
#endif