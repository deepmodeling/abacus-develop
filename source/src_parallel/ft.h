#ifndef FT_H
#define FT_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "parallel_pw.h"

#if defined __FFTW2
#include "fftw.h"
#elif defined __FFTW3
#include "fftw3.h"
#else
#include <fftw3-mpi.h>
//#include "fftw3-mpi_mkl.h"
#endif
typedef fftw_complex FFTW_COMPLEX;

class FFT: public Parallel_PW
{
public:

	FFT();
	~FFT();

	// mohan add 'const' 2021-02-25
	void FFT3D(std::complex<double> *psi,const int sign);
//	void FFT3D(double *psi, const int sign);
//	void FFT3D(matrix &psi, const int sign); // just for now

	const int& nx(void)const{return plan_nx;}
	const int& ny(void)const{return plan_ny;}
	const int& nz(void)const{return plan_nz;}
	const int& nxyz(void)const{return nxx;}
// #ifdef __MPI
	void setup_MPI_FFT3D(const int nx, const int ny, const int nz, const int nxx,const bool in_pool);
// #else
	void setupFFT3D(const int nx, const int ny, const int nz);
	void setupFFT3D_2(void);
// #endif
private:
	double scale_xyz;
	int plan_nx,plan_ny,plan_nz;
	int nxx;
	bool FFTWsetupwasdone;
	int test;
	std::complex<double> *aux4plan;
// #ifdef __MPI
	void P3DFFT(std::complex<double> *psi, const int sign);
	void fftxy(std::complex<double> *psi, const int sign);
	void fftz(std::complex<double> *psi_in, const int sign, std::complex<double> *psi_out);
	void scatter(std::complex<double> *psi, const int sign);
	std::complex<double> *aux;
	fftw_plan planplus_x;
	fftw_plan planplus_y; // mohan fix bug 2009-10-10
	fftw_plan planplus_z;
	fftw_plan planminus_x;
	fftw_plan planminus_y;
	fftw_plan planminus_z;
	double scale_xy;
	double scale_z;
	bool in_pool;
	int rank_use;
	int nproc_use;
	int *plane;
	int *sentc;
	int *recvc;
	int *sdis;
	int *rdis;
	int *sum;
// #else
	void SFFT3D(std::complex<double> *data, const int sign);

#if defined __FFTW2
	fftwnd_plan plus_plan;
	fftwnd_plan minus_plan;
#elif defined __FFTW3
	fftw_plan plus_plan;
	fftw_plan minus_plan;
#endif
// #endif
};

#endif
/*
 * This is a 3D FFT which uses the FFTW package.
 *
 * The transform performed is:
 *
 * if (sign == -1)
 * data(kx,ky,kz) <- sum_(x,y,z) { exp(-i*k.r))*data(x,y,z) } / (nx*ny*nz)
 *
 * if (sign == 1)
 * data(x,y,z) <- sum_(kx,ky,kz) { exp(i*k.r))*data(kx,ky,kz) }
 *
 * nx, ny, nz:  size of FFT box
 *
 * FFTW routines do not divide by (nx*ny*nz) in the sign==-1 case, so
 * I do it by hand.
 *
 * The plus_- and minus_plan's are static so that we do the setup only once
 * and on all future calls on FFT boxes of the same size we use those
 * plans.
 *
 * NOTE:  setupFFT3D() must be called before using FFT3D().
 *
 */

