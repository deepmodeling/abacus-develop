#include "fft_base.h"
#include <memory>
// #include "module_psi/psi.h"
#ifndef FFT_TEMP_H
#define FFT_TEMP_H
namespace ModulePW
{
class FFT_Bundle
{
    public:
        FFT_Bundle(){};
        ~FFT_Bundle();
        /**
         * @brief Constructor with device and precision.
         * @param device_in  device type, cpu or gpu.
         * @param precision_in  precision type, single or double.
         * 
         * the function will check the input device and precision, 
         * and set the device and precision.
         */
        FFT_Bundle(std::string device_in,std::string precision_in)
        :device(device_in),precision(precision_in){};
        
        /**
         * @brief Set device and precision.
         * @param device_in  device type, cpu or gpu.
         * @param precision_in  precision type, single or double.
         * 
         * the function will check the input device and precision, 
         * and set the device and precision.
         */
        void setfft(std::string device_in,std::string precision_in);

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

        void initfftmode(int fft_mode_in){this->fft_mode = fft_mode_in;}

        void setupFFT();

        void clearFFT();
        
        void clear();
        
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
        void fftzfor(std::complex<FPTYPE>* in, 
                     std::complex<FPTYPE>* out) const;
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
        void fftxyfor(std::complex<FPTYPE>* in, 
                     std::complex<FPTYPE>* out) const;
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
        void fftzbac(std::complex<FPTYPE>* in, 
                     std::complex<FPTYPE>* out) const;
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
        void fftxybac(std::complex<FPTYPE>* in, 
                     std::complex<FPTYPE>* out) const;
        
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
        void fftxyr2c(FPTYPE* in, 
                     std::complex<FPTYPE>* out) const;
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
        void fftxyc2r(std::complex<FPTYPE>* in, 
                     FPTYPE* out) const;

        template <typename FPTYPE, typename Device>
        void fft3D_forward(const Device* ctx, 
                          std::complex<FPTYPE>* in, 
                          std::complex<FPTYPE>* out) const;
        template <typename FPTYPE, typename Device>
        void fft3D_backward(const Device* ctx, 
                            std::complex<FPTYPE>* in, 
                            std::complex<FPTYPE>* out) const;

    private:
        int  fft_mode = 0; 
        bool float_flag=false;
        bool float_define=true;
        bool double_flag=false;
        std::shared_ptr<FFT_BASE<float>> fft_float=nullptr;
        std::shared_ptr<FFT_BASE<double>> fft_double=nullptr;
        
        std::string device = "cpu";
        std::string precision = "double";
};   
} // namespace ModulePW
#endif // FFT_H

