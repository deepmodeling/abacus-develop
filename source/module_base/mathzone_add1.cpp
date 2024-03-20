#include "mathzone_add1.h"
#include "timer.h"
#include "mathzone.h"
#include "global_variable.h"
#include "constants.h"
#include "global_function.h"
#include "math_sphbes.h"

#if defined __FFTW2
#include "fftw.h"
#elif defined __FFTW3
#include "fftw3.h" // mohan update 2021-05-06
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#else
#include <fftw3-mpi.h> 
//#include "fftw3-mpi_mkl.h"
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#endif
typedef fftw_complex FFTW_COMPLEX;

//#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleBase
{

double** Mathzone_Add1::c_ln_c = nullptr;
double** Mathzone_Add1::c_ln_s = nullptr;

Mathzone_Add1::Mathzone_Add1()
{}

Mathzone_Add1::~Mathzone_Add1()
{}


double Mathzone_Add1::factorial (const int& l)
{
	if (l == 0 || l == 1) return 1.0;	
	else return l*factorial(l-1);
}

double Mathzone_Add1::dualfac (const int& l)
{
	if (l == -1 || l == 0) return 1.0;
	else return l * dualfac (l-2); 
}


void Mathzone_Add1::SplineD2 // modified by pengfei 13-8-8 add second derivative as a  condition
(
 	const double *rad,
	const double *rad_f,
	const int& mesh,
	const double &yp1, // if yp1 > ypmax, consider the second derivative
	const double &ypn, 
	double* y2
)
{
	ModuleBase::timer::tick("Mathzone_Add1","SplineD2");
	
	double dx1, dx2, dy1, dy2, p, qn, sig, un;
	double *u;
	
	u = new double[mesh-1];
	const double ypmax = 99999.00;
	
	if (yp1 > ypmax) 
	{       
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else
	{
		y2[0] = -0.5;	
		
		dx1 = rad[1] - rad[0];
		dy1 = rad_f[1] - rad_f[0];

		u[0] = ( 3.0 / dx1 ) * (dy1 / dx1 - yp1);
	}
	
	for(int i = 1; i < mesh-1 ; i++)
	{
		dx1 = rad[i] - rad[i-1];
		dx2 = rad[i+1] - rad[i-1];
		dy1 = rad_f[i+1] - rad_f[i];
		dy2 = rad_f[i] - rad_f[i-1];
		
		sig = dx1 / dx2;

		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		
		u[i] = dy1 / (dx2 - dx1) - dy2 / dx1;
		u[i] = (6.0 * u[i] / dx2 - sig * u[i-1]) / p;
	}
	
	if (ypn > ypmax)
	{
		qn = un = 0.0;
	}
	else
	{
		dx1 = rad[mesh-1] - rad[mesh-2];
		dy1 = rad_f[mesh-1] - rad_f[mesh-2];
		
		qn = 0.5;
		un = 3.0 / dx1 * (ypn - dy1 / dx1);
	}

	y2[mesh-1] = (un - qn * u[mesh-2]) / (qn * y2[mesh-2] + 1.0);
	
	for(int i = mesh-2; i >= 0; i--)
	{
		y2[i] = y2[i] * y2[i+1] + u[i];
	}

	delete[] u;

	ModuleBase::timer::tick("Mathzone_Add1","SplineD2");
}

// Peize Lin add openmp 2019-12-13	
void Mathzone_Add1::Cubic_Spline_Interpolation
(
 	const double * const rad,
	const double * const rad_f,
	const double * const y2,
	const int& mesh,
	const double * const r,
	const int& rsize,
	double * const y,
	double * const dy
)
{	
	ModuleBase::timer::tick("Mathzone_Add1","Cubic_Spline_Interpolation");

#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for(int m = 0; m < rsize ; m++)
	{
		int klo = 0;
		int khi = mesh-1;
		while (khi - klo > 1)
		{
			const int k = (khi + klo) / 2 ;
			if(rad[k] > r[m]) khi = k;
			else klo = k;
		}

		const double h = rad[khi] - rad[klo];
		if(h == 0.0) 
		{
			std::cout << "Cubic_Spline_Interpolation, h == 0.0 so that cannot be divided" << std::endl;
			//ModuleBase::WARNING_QUIT("Cubic_Spline_Interpolation","h == 0.0 so that cannot be divided");
			exit(0);
		}

		const double a = (rad[khi] - r[m]) / h;
		const double b = (r[m] - rad[klo]) / h;
		
		const double dy_tmp = (rad_f[khi] - rad_f[klo]) / h - 
						(3.0 * a * a - 1.0) / 6.0 * h * y2[klo] + 
										( 3.0 * b * b - 1.0) / 6.0 * h * y2[khi];
		dy[m] = dy_tmp;
		const double y_tmp = a * rad_f[klo] + b * rad_f[khi] + ((a*a*a - a) * y2[klo] + (b*b*b - b) * y2[khi]) * (h*h) / 6.0;
		y[m] = y_tmp;
		//const double ddy_tmp = a * y2[klo] + b * y2 [khi];
		//ddy[m] = ddy_tmp;
	}

	ModuleBase::timer::tick("Mathzone_Add1","Cubic_Spline_Interpolation");
}

/// Interpolation for Numerical Orbitals
double Mathzone_Add1::Uni_RadialF
(
	const double* old_phi,
	const int& msh,
	const double& dr,
	const double& newr
)
{
  	double h1, h2, h3, f1, f2, f3, f4;
  	double g1, g2, x1, x2, y1, y2, f;
  	double result;
	double rmax = (msh-1) * dr;

  	if (newr < 0.0)
  	{  
		throw std::runtime_error("newr should >= 0. "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

  	}
  	else if ( rmax <= newr	) 
	{
		f = 0.0;
	}
  	else
  	{
		// the old position. 
		const int m = static_cast<int> (newr / dr) + 1;

       // Spline like interpolation
    	if (m == 1)
		{
      		h2 = dr;
      		h3 = dr;

      		f2 = old_phi[m-1];
      		f3 = old_phi[m];
      		f4 = old_phi[m+1];

      		h1 = -(h2+h3);
      		f1 = f4;
    	}
    	else if (m == (msh-1))
		{
      		h1 = dr;
      		h2 = dr;

      		f1 = old_phi[m-2];
      		f2 = old_phi[m-1];
      		f3 = old_phi[m];

      		h3 = -(h1+h2);
      		f4 = f1;
    	}
    	else
		{
			h1 = dr;
      		h2 = dr;
      		h3 = dr;

      		f1 = old_phi[m-2];
      		f2 = old_phi[m-1];
      		f3 = old_phi[m];
      		f4 = old_phi[m+1];
    	}

        // Calculate the value at newr

    	g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    	g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    	x1 = newr - (m-1)*dr;
   	 	x2 = newr - m*dr;
    	y1 = x1/h2;
    	y2 = x2/h2;

    	f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       	   + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  	} 
  	result = f;

  	return result;
}

void Mathzone_Add1::Uni_Deriv_Phi
(
	const double *radf,
	const int &mesh,
	const double &dr,
	const int &nd,
	double *phind
)
{
	ModuleBase::timer::tick("Mathzone_Add1", "Uni_Deriv_Phi");
	int FFT_NR = 2*mesh-1;  // FFT_NR = 16019
	// std::cout << "\n mesh=" << mesh << ", radf[8010]=" << radf[8010] <<  ", radf[8009]=" << radf[8009] ;
	// mesh=8010, radf[8010]=4.396478951532926e-01, radf[8009]=0.000000000000000e+00

    fftw_complex *fft_phir = new fftw_complex[FFT_NR];
    fftw_complex *fft_phik = new fftw_complex[FFT_NR];
    fftw_complex *fft_ndphik = new fftw_complex[FFT_NR];
    fftw_complex *fft_ndphir = new fftw_complex[FFT_NR];
	fftw_plan p1;
	fftw_plan p2;

	////CAREFUL: POINT 0 is OF GOOD IMPORTANCE	
	//for (int ir = 0; ir < FFT_NR/2; ++ir)
	//{
		//fft_phir[ir].re = radf[ir];
		//fft_phir[ir].im = 0.0;
	//}

	//for (int ir = FFT_NR/2; ir < FFT_NR; ++ir)
	//{
		//int jr = FFT_NR - ir;
		//fft_phir[ir].re = radf[jr];
		//fft_phir[ir].im = 0.0;
	//}
	//
	// second call: different value at [8010]; FFT_NR = 16019, FFT_NR/2=8009 :: 
	// CAREFUL: POINT 0 is OF GOOD IMPORTANCE	
	for (int ir = 0; ir < FFT_NR/2; ++ir)	   // ik = 0    1    ... 8008
	{
		c_re(fft_phir[ir]) = radf[ir];
		c_im(fft_phir[ir]) = 0.0;
	}

	for (int ir = FFT_NR/2; ir < FFT_NR; ++ir) // ir = 8009 8010 ... 16018
	{
		//int jr = FFT_NR - ir ;			   // jr = 8010 8009 ... 1  
		int jr = FFT_NR - ir -1 ;			//     ->  8009 8008 ... 0
		c_re(fft_phir[ir]) = radf[jr];
		c_im(fft_phir[ir]) = 0.0;
	}

	// FFTW
#if defined __FFTW3
	//std::cout << "\n Call FFTW3 ";
	p1 = fftw_plan_dft_1d(FFT_NR, fft_phir, fft_phik, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p1);
	//fftw_destroy_plan(p1);
#elif defined __FFTW2
	//std::cout << "\n Call FFTW2 ";
	p1 = fftw_create_plan(FFT_NR, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_one(p1, fft_phir, fft_phik);
	//fftw_destroy_plan(p1);
#endif

	
	double dk_uniform = ModuleBase::TWO_PI / FFT_NR / dr;
	
	//for (int ik = 0; ik < FFT_NR/2; ik++)
	//{
		//double kp = ik * dk_uniform;
		//fft_ndphik[ik].re = pow(kp, nd) * fft_phik[ik].re;
		//fft_ndphik[ik].im = 0.0;
	//}
	
	//for (int ik = FFT_NR/2; ik < FFT_NR; ik++)
	//{
		//double kp = -(FFT_NR - ik)* dk_uniform;
		//fft_ndphik[ik].re = pow(kp, nd) * fft_phik[ik].re;
		//fft_ndphik[ik].im = 0.0;
	//}

	for (int ik = 0; ik < FFT_NR/2; ik++)				// ik = 0    1    ... 8008
	{
		double kp = ik * dk_uniform;
		c_re(fft_ndphik[ik]) = pow(kp, nd) * c_re(fft_phik[ik]);
		c_im(fft_ndphik[ik]) = 0.0;
	}
	for (int ik = FFT_NR/2; ik < FFT_NR; ik++)			// ik = 8009 8010 ... 16018
	{
		//double kp = -(FFT_NR - ik )* dk_uniform;		//(...)  = 8010 8009 ... 1  
		double kp = -(FFT_NR - ik -1)* dk_uniform;		//(...) -> 8009 8008 ... 0
		c_re(fft_ndphik[ik]) = pow(kp, nd) * c_re(fft_phik[ik]);
		c_im(fft_ndphik[ik]) = 0.0;
	}

#if defined __FFTW3
	p2 = fftw_plan_dft_1d(FFT_NR, fft_ndphik, fft_ndphir, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	//fftw_destroy_plan(p2);
#elif defined __FFTW2
	p2 = fftw_create_plan(FFT_NR, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_one(p2, fft_ndphik, fft_ndphir);
	//fftw_destroy_plan(p2);
#endif

	bool is_re;
	double fac;
	if (nd % 4 == 0) 
	{
		is_re = true; 
		fac = 1.0;
	}
	else if (nd % 4 == 1)
	{
		is_re = false;
		fac = -1.0;
	}
	else if (nd % 4 == 2)
	{
		is_re = true;
		fac = -1.0;
	}
	else
	{
		is_re = false;
		fac = 1.0;
	}

	for (int ir = 0; ir < mesh; ir++)
	{
		if (is_re)
		{
			phind[ir] = fac * c_re(fft_ndphir[ir]) / FFT_NR;
		}
		else
		{
			phind[ir] = fac * c_im(fft_ndphir[ir]) / FFT_NR;
		}
	}
	
	fftw_destroy_plan (p1);
	fftw_destroy_plan (p2);

    delete [] fft_phir;
    delete [] fft_phik;
    delete [] fft_ndphik;
    delete [] fft_ndphir;
	
	ModuleBase::timer::tick("Mathzone_Add1", "Uni_Deriv_Phi");
}

}
