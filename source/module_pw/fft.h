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
	void initfft(int nx_in, int bigny_in, int nz_in, int liy_in, int riy_in, int ns_in, int nplane_in, 
				 int nproc_in, bool gamma_only_in, bool mpifft_in = false);
	void setupFFT();
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
	void initplan();
	void initplan_mpi();
#ifdef __MIX_PRECISION
	void initplanf();
	void initplanf_mpi();
#endif
	
public:
	int nx,ny,nz;
	int nxy;
	int bigny;
	int bignxy;
	int liy,riy;// liy: the left edge of the pw ball in the y direction; riy: the right edge of the pw ball in the y direction
	int ns; //number of sticks
	int nplane; //number of x-y planes
	int maxgrids; // max between nz * ns and bignxy * nplane
	int nproc; // number of proc.
	std::complex<double> *aux1, *aux2; //fft space, [maxgrids]
	double *r_rspace; //real number space for r, [nplane * nx *ny]
#ifdef __MIX_PRECISION
	std::complex<float> *auxf1, *auxf2; //fft space, [maxgrids]
	float *rf_rspace; //real number space for r, [nplane * nx *ny]
#endif


private:
	bool gamma_only;
	bool destroyp;
	bool mpifft; // if use mpi fft, only used when define __FFTW3_MPI
	// fftw_plan plan2r2c;
	// fftw_plan plan2c2r;
	// fftw_plan plan1for;
	// fftw_plan plan1bac;
	// fftw_plan plan2for;
	// fftw_plan plan2bac;
	fftw_plan planzfor;
	fftw_plan planzbac;
	fftw_plan planxfor1;
	fftw_plan planxbac1;
	fftw_plan planxfor2;
	fftw_plan planxbac2;
	fftw_plan planyfor;
	fftw_plan planybac;
	fftw_plan planyr2c;
	fftw_plan planyc2r;
#ifdef __MIX_PRECISION
	bool destroypf;
	fftwf_plan planfzfor;
	fftwf_plan planfzbac;
	fftwf_plan planfxfor1;
	fftwf_plan planfxbac1;
	fftwf_plan planfxfor2;
	fftwf_plan planfxbac2;
	fftwf_plan planfyfor;
	fftwf_plan planfybac;
	fftwf_plan planfyr2c;
	fftwf_plan planfyc2r;
#endif

};
}

#endif

