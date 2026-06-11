#ifndef FFT_TEMP_H
#define FFT_TEMP_H

#include "fft_base.h"
#include "fft_cpu.h"

#include <memory>
namespace ModuleBase
{
class FFT_Bundle
{
  public:
    FFT_Bundle() {};
    ~FFT_Bundle();
    /**
     * @brief Constructor with device and precision.
     * @param device_in  device type, cpu or gpu.
     * @param precision_in  precision type, single or double.
     *
     * the function will check the input device and precision,
     * and set the device and precision.
     */
    FFT_Bundle(std::string device_in, std::string precision_in) : device(device_in), precision(precision_in) {};

    /**
     * @brief Set device and precision.
     * @param device_in  device type, cpu or gpu.
     * @param precision_in  precision type, single or double.
     *
     * the function will check the input device and precision,
     * and set the device and precision.
     */
    void setfft(std::string device_in, std::string precision_in);

    /**
     * @brief Set the DSP cluster id for the FFT_DSP backend.
     * @param id  cluster id, typically computed as (MPI rank % dsp_count).
     *
     * Caller-injected DSP routing info; only used when device == "dsp".
     */
    void set_dsp_cluster_id(int id) { this->dsp_cluster_id_ = id; }

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
     *
     * the function will initialize the many-fft parameters
     * Wheatley in cpu or gpu device.
     */
    void initfft(int nx_in,
                 int ny_in,
                 int nz_in,
                 int lixy_in,
                 int rixy_in,
                 int ns_in,
                 int nplane_in,
                 int nproc_in,
                 bool gamma_only_in,
                 bool xprime_in = true,
                 bool mpifft_in = false);

    /**
     * @brief Initialize the fft mode.
     * @param fft_mode_in  fft mode.
     *
     * the function will initialize the fft mode.
     */

    void initfftmode(int fft_mode_in)
    {
        this->fft_mode = fft_mode_in;
    }

    /**
     * @brief Initialize the batch FFT size.
     * @param batch_size_in  batch size for batch FFT (1-128)
     *
     * the function will initialize the batch FFT size.
     */
    void init_batch_size(int batch_size_in)
    {
        this->batch_size = batch_size_in;
    }

    void setupFFT();

    void clearFFT();

    void clear();

    void resource_handler(const int flag) const;
    /**
     * @brief Get the real space data.
     * @return FPTYPE*  the real space data.
     *
     * the function will return the real space data,
     * which is used in the cpu-like fft.
     */
    template <typename FPTYPE>
    FPTYPE* get_rspace_data() const;
    /**
     * @brief Get the auxr data.
     * @return std::complex<FPTYPE>*  the auxr data.
     *
     * the function will return the auxr data,
     * which is used in the cpu-like fft.
     */
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxr_data() const;
    /**
     * @brief Get the auxg data.
     * @return std::complex<FPTYPE>*  the auxg data.
     *
     * the function will return the auxg data,
     * which is used in the cpu-like fft.
     */
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxg_data() const;
    /**
     * @brief Get the auxr 3d data.
     * @return std::complex<FPTYPE>*  the auxr 3d data.
     *
     * the function will return the auxr 3d data,
     * which is used in the gpu-like fft.
     */
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxr_3d_data() const;

    /**
     * @brief Forward fft in z direction.
     * @param in  input data.
     * @param out  output data.
     *
     * The function will do the forward many fft in z direction,
     * As an interface, the function will call the fftzfor in the
     * accurate fft class.
     * which is used in the cpu-like fft.
     */
    template <typename FPTYPE>
    void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    /**
     * @brief Forward fft in x-y direction.
     * @param in  input data.
     * @param out  output data.
     *
     * the function will do the forward fft in x and y direction,
     * which is used in the cpu-like fft.As an interface,
     * the function will call the fftxyfor in the accurate fft class.
     */
    template <typename FPTYPE>
    void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    /**
     * @brief Backward fft in z direction.
     * @param in  input data.
     * @param out  output data.
     *
     * the function will do the backward many fft in z direction,
     * which is used in the cpu-like fft.As an interface,
     * the function will call the fftzbac in the accurate fft class.
     */
    template <typename FPTYPE>
    void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    /**
     * @brief Backward fft in x-y direction.
     * @param in  input data.
     * @param out  output data.
     *
     * the function will do the backward fft in x and y direction,
     * which is used in the cpu-like fft.As an interface,
     * the function will call the fftxybac in the accurate fft class.
     */
    template <typename FPTYPE>
    void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    /**
     * @brief Real to complex fft in x-y direction.
     * @param in  input data.
     * @param out  output data.
     *
     * the function will do the real to complex fft in x and y direction,
     * which is used in the cpu-like fft.As an interface,
     * the function will call the fftxyr2c in the accurate fft class.
     */
    template <typename FPTYPE>
    void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const;
    /**
     * @brief Complex to real fft in x-y direction.
     * @param in  input data.
     * @param out  output data.
     *
     * the function will do the complex to real fft in x and y direction,
     * which is used in the cpu-like fft.As an interface,
     * the function will call the fftxyc2r in the accurate fft class.
     */
    template <typename FPTYPE>
    void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const;

    template <typename FPTYPE>
    void fft3D_forward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fft3D_backward(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

    // Batch FFT methods
    /**
     * @brief Setup batch FFT plans and allocate batch buffers
     *
     * Initializes cufftPlanMany for batch FFT operations on GPU.
     * No-op on CPU devices. Batch execution APIs throw if batch FFT is unavailable.
     */
    void setupBatchFFT();

    /**
     * @brief Forward batch 3D FFT
     * @param in_batch  input data batch
     * @param out_batch output data batch
     * @param batch_count actual number of FFTs to process
     *
     * Performs batch_count forward 3D FFTs in a single kernel launch (GPU).
     * Falls back to sequential FFTs on CPU or if batch FFT not available.
     */
    template <typename FPTYPE>
    void fft3D_forward_batch(std::complex<FPTYPE>* in_batch,
                             std::complex<FPTYPE>* out_batch,
                             int batch_count) const;

    /**
     * @brief Backward batch 3D FFT
     * @param in_batch  input data batch
     * @param out_batch output data batch
     * @param batch_count actual number of FFTs to process
     *
     * Performs batch_count backward 3D FFTs in a single kernel launch (GPU).
     * Falls back to sequential FFTs on CPU or if batch FFT not available.
     */
    template <typename FPTYPE>
    void fft3D_backward_batch(std::complex<FPTYPE>* in_batch,
                              std::complex<FPTYPE>* out_batch,
                              int batch_count) const;

    /**
     * @brief Check if batch FFT is available for given precision
     * @return true if batch FFT is ready and available
     */
    template <typename FPTYPE>
    bool is_batch_fft_available() const;

    /**
     * @brief Get maximum batch size
     * @return maximum batch size for batch FFT operations
     */
    template <typename FPTYPE>
    int get_batch_size() const;

    /**
     * @brief Get batch input buffer
     * @return pointer to device batch input buffer
     */
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_batch_input_buffer() const;

    /**
     * @brief Get batch output buffer
     * @return pointer to device batch output buffer
     */
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_batch_output_buffer() const;

  private:
    int fft_mode = 0;
    int batch_size = 8;  // Default batch size for batch FFT
    bool float_flag = false;
    bool double_flag = false;

    // Primary FFT objects (CPU or GPU depending on device setting)
    std::shared_ptr<FFT_BASE<float>> fft_float = nullptr;
    std::shared_ptr<FFT_BASE<double>> fft_double = nullptr;

    // CPU FFT objects for fallback when device="gpu"
    // These are used by non-templated CPU-style FFT operations (get_auxg_data, fftxyfor, etc.)
    std::shared_ptr<FFT_BASE<float>> fft_float_cpu = nullptr;
    std::shared_ptr<FFT_BASE<double>> fft_double_cpu = nullptr;

    std::string device = "cpu";
    std::string precision = "double";
    int dsp_cluster_id_ = 0;
};
// Use RAII (Resource Acquisition Is Initialization) to 
// control the resources used by hthread when setting the DSP
struct FFT_Guard
  {
      const FFT_Bundle& fft_;
      FFT_Guard(const FFT_Bundle& fft) : fft_(fft) 
        {fft_.resource_handler(1);}
      ~FFT_Guard()
      {
        fft_.resource_handler(0);
      }
  };

} // namespace ModuleBase
#endif // FFT_H
