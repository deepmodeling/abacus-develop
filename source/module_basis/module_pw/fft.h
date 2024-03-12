#ifndef FFT_H
#define FFT_H

#include <complex>
#include <string>

#include "fftw3.h"
#if defined(__FFTW3_MPI) && defined(__MPI)
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif

#if defined(__CUDA) || defined(__UT_USE_CUDA)
#include "cufft.h"
#include "cuda_runtime.h"
#endif

#if defined(__ROCM) || defined(__UT_USE_ROCM)
#include <hipfft/hipfft.h>
#include <hip/hip_runtime.h>
#endif

//Temporary: we donot need psi. However some GPU ops are defined in psi, which should be moved into module_base or module_gpu
#include "module_psi/psi.h"
// #ifdef __ENABLE_FLOAT_FFTW
// #include "fftw3f.h"
// #if defined(__FFTW3_MPI) && defined(__MPI)
// #include "fftw3f-mpi.h"
// //#include "fftw3-mpi_mkl.h"
// #endif
// #endif

namespace ModulePW
{

class FFT
{
public:

	FFT();
	~FFT();
	void clear(); //reset fft
	
	// init parameters of fft
	void initfft(int nx_in, int ny_in, int nz_in, int lixy_in, int rixy_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool xprime_in = true, bool mpifft_in = false);

	//init fftw_plans
	void setupFFT(); 

	//destroy fftw_plans
	void cleanFFT();

#if defined(__ENABLE_FLOAT_FFTW)
	void cleanfFFT();
#endif // defined(__ENABLE_FLOAT_FFTW)

    template <typename FPTYPE>
    void fftzfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fftzbac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fftxyfor(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fftxybac(std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fftxyr2c(FPTYPE* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE>
    void fftxyc2r(std::complex<FPTYPE>* in, FPTYPE* out) const;

    template <typename FPTYPE, typename Device>
    void fft3D_forward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;
    template <typename FPTYPE, typename Device>
    void fft3D_backward(const Device* ctx, std::complex<FPTYPE>* in, std::complex<FPTYPE>* out) const;

  public:
	//init fftw_plans
	void initplan(const unsigned int& flag = 0);
	// We have not support mpi fftw yet.
	// void initplan_mpi();
	//init fftwf_plans
#if defined(__ENABLE_FLOAT_FFTW)
	void initplanf(const unsigned int& flag = 0);
#endif // defined(__ENABLE_FLOAT_FFTW)
	// void initplanf_mpi();

public:
	int fftnx=0, fftny=0;
	int fftnxy=0;
	int ny=0, nx=0, nz=0;
	int nxy=0;
	bool xprime = true; // true: when do recip2real, x-fft will be done last and when doing real2recip, x-fft will be done first; false: y-fft
                         // For gamma_only, true: we use half x; false: we use half y
	int lixy=0,rixy=0;// lixy: the left edge of the pw ball in the y direction; rixy: the right edge of the pw ball in the x or y direction
	int ns=0; //number of sticks
	int nplane=0; //number of x-y planes
	int nproc=1; // number of proc.

    template <typename FPTYPE>
    FPTYPE* get_rspace_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxr_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxg_data() const;
    template <typename FPTYPE>
    std::complex<FPTYPE>* get_auxr_3d_data() const;

	int fft_mode = 0; ///< fftw mode 0: estimate, 1: measure, 2: patient, 3: exhaustive 

  private:
    bool gamma_only = false;
    bool mpifft = false; // if use mpi fft, only used when define __FFTW3_MPI

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
//	fftw_plan plan3dforward;
//	fftw_plan plan3dbackward;

#if defined(__CUDA)
    cufftHandle c_handle = {};
    cufftHandle z_handle = {};
#elif defined(__ROCM)
    hipfftHandle c_handle = {};
    hipfftHandle z_handle = {};
#endif

#if defined(__ENABLE_FLOAT_FFTW)
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
#endif // defined(__ENABLE_FLOAT_FFTW)

    mutable std::complex<float>* c_auxr_3d = nullptr;  // fft space
    mutable std::complex<double>* z_auxr_3d = nullptr; // fft space

    mutable std::complex<float>*c_auxg = nullptr, *c_auxr = nullptr;  // fft space,
    mutable std::complex<double>*z_auxg = nullptr, *z_auxr = nullptr; // fft space

    mutable float* s_rspace = nullptr;  // real number space for r, [nplane * nx *ny]
    mutable double* d_rspace = nullptr; // real number space for r, [nplane * nx *ny]

    std::string device = "cpu";
    std::string precision = "double";

public:
    void set_device(std::string device_);
    void set_precision(std::string precision_);

};
}

#endif

