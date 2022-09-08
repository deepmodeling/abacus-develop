#ifndef FFT_H
#define FFT_H

#include <complex>
#include <string>

#include "fftw3.h"
#if defined(__FFTW3_MPI) && defined(__MPI)
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif

// #ifdef __MIX_PRECISION
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

	void fftzfor(std::complex<double>* & in, std::complex<double>* & out);
	void fftzbac(std::complex<double>* & in, std::complex<double>* & out);
	void fftxyfor(std::complex<double>* & in, std::complex<double>* & out);
	void fftxybac(std::complex<double>* & in, std::complex<double>* & out);
	void fftxyr2c(double * &in, std::complex<double>* & out);
	void fftxyc2r(std::complex<double>* & in, double* & out);

#ifdef __MIX_PRECISION
	void cleanfFFT();
	void fftfzfor(std::complex<float>* & in, std::complex<float>* & out);
	void fftfzbac(std::complex<float>* & in, std::complex<float>* & out);
	void fftfxyfor(std::complex<float>* & in, std::complex<float>* & out);
	void fftfxybac(std::complex<float>* & in, std::complex<float>* & out);
	void fftfxyr2c(float * &in, std::complex<float>* & out);
	void fftfxyc2r(std::complex<float>* & in, float* & out);
#endif

public:
	//init fftw_plans
	void initplan(); 
	void initplan_mpi();
#ifdef __MIX_PRECISION
	//init fftwf_plans
	void initplanf(); 
	void initplanf_mpi();
#endif
	
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
	std::complex<double> *auxg=nullptr, *auxr=nullptr; //fft space
	double *r_rspace=nullptr; //real number space for r, [nplane * nx *ny]
#ifdef __MIX_PRECISION
	std::complex<float> *auxfg=nullptr, *auxfr=nullptr; //fft space,
	float *rf_rspace=nullptr; //real number space for r, [nplane * nx *ny]
#endif


private:
	bool gamma_only=false;
	bool mpifft=false; // if use mpi fft, only used when define __FFTW3_MPI

	bool destroyp=true;
	fftw_plan planzfor;
	fftw_plan planzbac;
	fftw_plan planxfor1;
	fftw_plan planxbac1;
	fftw_plan planxfor2;
	fftw_plan planxbac2;
	fftw_plan planyfor;
	fftw_plan planybac;
	fftw_plan planxr2c;
	fftw_plan planxc2r;
	fftw_plan planyr2c;
	fftw_plan planyc2r;
#ifdef __MIX_PRECISION
	bool destroypf=true;
	fftwf_plan planfzfor;
	fftwf_plan planfzbac;
	fftwf_plan planfxfor1;
	fftwf_plan planfxbac1;
	fftwf_plan planfxfor2;
	fftwf_plan planfxbac2;
	fftwf_plan planfyfor;
	fftwf_plan planfybac;
	fftwf_plan planfxr2c;
	fftwf_plan planfxc2r;
	fftwf_plan planfyr2c;
	fftwf_plan planfyc2r;
#endif

};
}

#endif

