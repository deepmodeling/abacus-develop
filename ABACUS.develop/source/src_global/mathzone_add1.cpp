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
#include <omp.h>

bool Mathzone_Add1::flag_jlx_expand_coef = false;
double** Mathzone_Add1::c_ln_c = nullptr;
double** Mathzone_Add1::c_ln_s = nullptr;

Mathzone_Add1::Mathzone_Add1()
{}

Mathzone_Add1::~Mathzone_Add1()
{}


/**********************************************
 * coefficients to expand jlx using
 * cos and sin, from SIESTA
 * *******************************************/
void Mathzone_Add1::expand_coef_jlx()
{
	timer::tick("Mathzone_Add1","expand_coef_jlx");

	int ir, il, in;
	const int L = 20;
	
	//allocation
	delete[] c_ln_c;
	delete[] c_ln_s;

	c_ln_c = new double*[L+1];
	c_ln_s = new double*[L+1];

	for(ir = 0; ir < L+1; ir++)
	{
		c_ln_c[ir] = new double[ir+1];
		c_ln_s[ir] = new double[ir+1];
	}

	//calculate initial value
	c_ln_c[0][0] = 0.0;
	c_ln_s[0][0] = 1.0;

	c_ln_c[1][0] = 0.0;
	c_ln_c[1][1] = -1.0;
	c_ln_s[1][0] = 1.0;
	c_ln_s[1][1] = 0.0;

	//recursive equation
	for(il = 2; il < Mathzone_Add1::sph_lmax+1; il++)
	{
		for(in = 0; in < il + 1; in++)
		{
			if(in  >= 2)
			{
				if(in == il)
				{
					c_ln_c[il][in] = -c_ln_c[il-2][in-2];
					c_ln_s[il][in] = -c_ln_s[il-2][in-2];
				}
				else
				{
					c_ln_c[il][in] = (2*il-1)*c_ln_c[il-1][in] - c_ln_c[il-2][in-2];
					c_ln_s[il][in] = (2*il-1)*c_ln_s[il-1][in] - c_ln_s[il-2][in-2];
				}
			}
			else
			{
				//in = 0 or 1, but here il = 2, so in != il
				c_ln_c[il][in] = (2*il-1)*c_ln_c[il-1][in];
				c_ln_s[il][in] = (2*il-1)*c_ln_s[il-1][in];
			}
		}
	}

	timer::tick("Mathzone_Add1","expand_coef_jlx");
	return;
}

void Mathzone_Add1::Spherical_Bessel
(           
	const int &msh, //number of grid points
	const double *r,//radial grid
	const double &q,    //
	const int &l,   //angular momentum
	double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
	double *sjp
)
{
	timer::tick ("Mathzone_Add1","Spherical_Bessel");
	
	assert (l <= Mathzone_Add1::sph_lmax);
	
	//creat coefficients
	if (!flag_jlx_expand_coef)
	{
		Mathzone_Add1::expand_coef_jlx ();
		flag_jlx_expand_coef = true;
	}
	
	//epsilon
	const double eps = 1.0E-10;
	
	/******************************************************************
	jlx = \sum_{n=0}^{l}(c_{ln}^{s}sin(x) + c_{ln}^{c}cos(x))/x^{l-n+1}
	jldx = \sum_{n=0}^{l}(c_{ln}^{s}sin(x) + c_{ln}^{c}cos(x))/x^{l-n+1}
	*******************************************************************/
	for (int ir = 0; ir < msh; ir++)
	{
		double qr = q * r[ir];
	
		//judge qr approx 0
		if (fabs(qr) < eps) 
		{
			if (l == 0) sj[ir] = 1.0;
			else sj[ir] = 0.0;

			if (l == 1) sjp[ir] = 1.0 / 3.0;
			else sjp[ir] = 0.0;
		}
		else
		{
			sj[ir] = 0.0;
			sjp[ir] = 0.0;
		
			double lqr = pow (qr, l);
			
			//divided fac
			double xqr = 1.0;
			
			for (int n = 0; n <= l; n++)
			{
				double com1 = (c_ln_s[l][n] * sin(qr)  + c_ln_c[l][n] * cos(qr)) * xqr;
				double com2 = (c_ln_s[l][n] * cos(qr) - c_ln_c[l][n] * sin(qr)) * xqr * qr
							+	(c_ln_s[l][n] * sin(qr) + c_ln_c[l][n] * cos(qr)) * n * xqr;
				
				sj[ir] += com1;
				sjp[ir] += (com2 - l*com1); 
			
				xqr *= qr;
			}

			sj[ir] /= lqr;
			sj[ir] /= (lqr * qr);
		}
	}

	timer::tick ("Mathzone_Add1","Spherical_Bessel");
	
	return;

}

double Mathzone_Add1::uni_simpson
(
 	const double* func,
	const int& mshr,
	const double& dr
)
{
	timer::tick("Mathzone_Add1","uni_simpson");
	
	assert(mshr >= 3);
	
	int ir=0;
	int idx=0;
	int msh_left=0;
	double sum = 0.0;
	
	//f(a)
	sum += func[0];
	
	//simpson 3/8 rule
	if(mshr % 2 == 0) 
	{
		//simpson 3/8 rule
		sum += 9.0 / 8.0 * (func[mshr-1] + 3.0*(func[mshr-2] + func[mshr-3]) + func[mshr-4]);
		msh_left = mshr - 3;
	}
	else 
	{
		msh_left = mshr;
	}
	
	//f(b)
	sum += func[msh_left-1];
	
	//points left
	for(ir = 0; ir < msh_left / 2; ir++)
	{
		idx = 2*ir+1;
		sum += 4.0 * func[idx];
	}
	
	for(ir = 1; ir < msh_left / 2; ir++)
	{
		idx = 2*ir;
		sum += 2.0 * func[idx];
	}

	timer::tick("Mathzone_Add1","uni_simpson");
	return sum * dr / 3.0;
}

void Mathzone_Add1::uni_radfft
(
	const int& n,
	const int& aml,
	const int& mshk,
	const double* arr_k,
	const int& mshr,
	const double* ri,
	const double& dr,
	const double* phir,
	double* phik
)
{
	timer::tick ("Mathzone_Add1","uni_radfft");
	
	//allocate memory
	double* fi = new double[mshr];
	double* fi_cp = new double[mshr];
	double* jl = new double[mshr];
	
	//function to be integrated: r^2 phir
	for (int ir = 0; ir < mshr; ir++)
	{
		fi_cp[ir] = pow(ri[ir], n) * phir[ir];
	}

	//integration
	for (int ik = 0; ik < mshk; ik++)
	{
		//calculate spherical bessel
		Sphbes::Spherical_Bessel(mshr, ri, arr_k[ik], aml, jl);
		
		//functions to be integrated
		for (int ir = 0; ir < mshr; ir++)
		{
			fi[ir] = fi_cp[ir] * jl[ir];
		}

		phik[ik] = uni_simpson (fi, mshr, dr);
	}

	//deallocate memory
	delete[] fi;		
	delete[] fi_cp;
	delete[] jl;

	timer::tick("Mathzone_Add1","uni_radfft");

	return;
}

void Mathzone_Add1::Sph_Bes (double x, int lmax, double *sb, double *dsb) 
{
	timer::tick("Mathzone_Add1","Spherical_Bessel");
	
  	int m, n, nmax;
  	double j0, j1, sf, tmp, si, co, ix;	// j0p, j1p, ix2

  	if (x < 0.0)
	{
		cout << "\nminus x is invalid for Spherical_Bessel" << endl;  
		exit(0); // mohan add 2021-05-06
	}
	 
  	const double xmin = 1E-10;

  	/* find an appropriate nmax */
  	nmax = lmax + 3*static_cast<int>(x) + 20;
  	if (nmax < 100) nmax = 100; 

  	/* allocate tsb */
  	double* tsb = new double[nmax+1];
  
  	/* if x is larger than xmin */

  	if (xmin < x)
	{
    	/* initial values*/
		tsb[nmax]   = 0.0;
    	tsb[nmax-1] = 1.0e-14;

    	/* downward recurrence from nmax-2 to lmax+2 */

    	for (n = nmax-1; (lmax+2) < n; n--)
		{
      		tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      		if (1.0e+250 < tsb[n-1])
			{
        		tmp = tsb[n-1];        
        		tsb[n-1] /= tmp;
        		tsb[n  ] /= tmp;
      		}
    	}

    	/* downward recurrence from lmax+1 to 0 */
		n = lmax + 3;
    	tmp = tsb[n-1];        
    	tsb[n-1] /= tmp;
    	tsb[n  ] /= tmp;

    	for (n = lmax+2; 0 < n; n--)
		{
      		tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      		if (1.0e+250 < tsb[n-1])
			{
        		tmp = tsb[n-1];
        		for (m = n-1; m <= lmax+1; m++)
				{
          			tsb[m] /= tmp;
        		}
      		}
    	}

    	/* normalization */
    	si = sin(x);
  	  	co = cos(x);
    	ix = 1.0/x;
    	//ix2 = ix*ix;
    	j0 = si*ix;
    	j1 = si*ix*ix - co*ix;

    	if (fabs(tsb[1]) < fabs(tsb[0])) sf = j0/tsb[0];
    	else                           sf = j1/tsb[1];

    	/* tsb to sb */
//	   	for (n = 0; n <= lmax+1; n++)
		for (n = 0; n <= lmax; n++)
		{
      		sb[n] = tsb[n]*sf;
    	}
		
    	/* derivative of sb */
    	dsb[0] = co*ix - si*ix*ix;
 //   	for (n = 1; n <= lmax; n++)
		for (n = 1; n < lmax; n++)
		{
      		dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    	}

		n = lmax;
      	dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sf*tsb[n+1] )/(2.0*(double)n + 1.0);
	} 

  	/* if x is smaller than xmin */
  	else 
	{
    	/* sb */
		for (n = 0; n <= lmax; n++ )
		{
      		sb[n] = 0.0;
    	}
    	sb[0] = 1.0;

    	/* derivative of sb */
    	dsb[0] = 0.0;
		dsb[1] = 1.0 / 3.0;
    	for (n = 2; n <= lmax; n++)
		{
//      	dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    		dsb[n] = 0.0;
		}
  	}

  	/* free tsb */
  	delete [] tsb;

  	timer::tick("Mathzone_Add1","Spherical_Bessel");

  	return;
}

void Mathzone_Add1::Sbt_new 
(
	const int& polint_order,
	const int& l,
	const double* k,
	const double& dk,
	const int& mshk,
	const double* r,
	const double& dr,
	const int& mshr,
	const double* fr,
	const int& rpow,
	double* fk
)
{
	timer::tick ("Mathzone_Add1","Sbt_new");

	//check parameter
	assert (mshr >= 1);
	assert (l >= 0);
	assert (rpow >= 0 && rpow <=2);
	
	//step 0
	//l is odd or even
	bool parity_flag;
	if (l % 2 == 0) 
	{
		parity_flag = true;
	}
	else 
	{
		parity_flag = false;
	}
	
	ZEROS (fk, mshk);
	
	if (polint_order != 3)
	{
		cout << "\nhigh order interpolation is not available!" << endl;
		exit(0); // mohan add 2021-05-06
		//QUIT();
	}
	
	/**********************************
	function multiplied by power of r
	for different polint_order
	**********************************/
	double* fr2;
	double* fr3;

	//polint_order == 1
	fr2 = new double[mshr];
	if (rpow == 0) for (int ir = 0; ir < mshr; ir++) fr2[ir] = fr[ir]*r[ir]*r[ir];
	else if (rpow == 1)	for (int ir = 0; ir < mshr; ir++) fr2[ir] = fr[ir]*r[ir];
	else if (rpow == 2) for (int ir = 0; ir < mshr; ir++) fr2[ir] = fr[ir];
	
	fr3 = new double[mshr];
	for (int ir = 0; ir < mshr; ir++) fr3[ir] = fr2[ir] * r[ir];

//	const int polint_order = 3;
	int nu_pol_coef  = (polint_order+1)*(mshk-1);
	double* polint_coef = new double[nu_pol_coef]; 

	//step 1
	//start calc
	if (parity_flag)
	{
		//even
		const int n = l/2;
		
		//coef for interpolation
		int ct = 0;

		double ft_save = fourier_cosine_transform (fr2, r, mshr, dr, k[0]);
		double dft_save = -fourier_sine_transform (fr3, r, mshr, dr, k[0]);
		
		for (int ik = 0; ik < mshk-1; ik++)
		{
			double ft0 = ft_save;
			double dft0 = dft_save;
			
			double ft1 = fourier_cosine_transform (fr2, r, mshr, dr, k[ik+1]);
			double dft1 = -fourier_sine_transform (fr3, r, mshr, dr, k[ik+1]);

			//double d2k = dk*dk;
			//double d3k = d2k*dk;
				
			double c0 = ft0;
			double c1 = dft0;
			double c2 = 3.0*(ft1-ft0)/dk/dk-(dft1+2.0*dft0)/dk;
			double c3 = (-2.0*(ft1-ft0)/dk+(dft1+dft0))/dk/dk;
				
			double k2 = k[ik]*k[ik];
			double k3 = k2*k[ik];
					
			polint_coef[ct] = c0-c1*k[ik]+c2*k2-c3*k3;
			polint_coef[ct+1] = c1-2.0*c2*k[ik]+3.0*c3*k2;
			polint_coef[ct+2] = c2-3.0*k[ik]*c3;
			polint_coef[ct+3] = c3;
			
			//test
	//		double x = (k[ik]+k[ik+1])/2;
	//		double x = k[ik];
//			double tmp = polint_coef[ct]+x*polint_coef[ct+1]+x*x*polint_coef[ct+2]+x*x*x*polint_coef[ct+3];
	//		double tmp_ana = fourier_cosine_transform (fr2, r, mshr, dr, x);
	//		cout << "\ninterp = " << tmp << " ana = " << tmp_ana << " diff = " << log(fabs(tmp-tmp_ana))/log(10);
			
			//update
			ct += (polint_order+1);
			ft_save = ft1;
			dft_save = dft1;
		}
		
		//store coefficients for calculation
		double* coef = new double[n+1];
		double fac = dualfac (l-1) / dualfac (l);

		for (int j = 0; j < n; j++)
		{
			coef[j] = fac;
			
			//update
			int twoj = 2*j;
			fac *= -static_cast<double>((l+twoj+1)*(l-twoj))/(twoj+2)/(twoj+1);
		}
		coef[n] = pow (-1.0, n) * dualfac (2*l-1) / factorial (l);

		//start calc
		//special case k = 0;
		if (n ==0) fk[0] = uni_simpson (fr2, mshr, dr); 
		else fk[0] = 0.0;

		//k > 0
		for (int j =0; j <=n; j++)
		{
			//initialize Snm
			double Snm = 0.0;
			for (int ik = 1; ik < mshk; ik++)
			{
				Snm += pol_seg_int (polint_order, polint_coef, 2*j, k, ik);
				double k2j = pow (k[ik], 2*j+1);

				fk[ik] += coef[j] * Snm / k2j;
			}
		}
		//free
		delete [] coef;
	}
	else
	{
		//odd
		const int n = (l-1)/2;
		
		//coef for interpolation
		int ct = 0;
		double ft_save, dft_save;
		ft_save = fourier_sine_transform (fr2, r, mshr, dr, k[0]);
		dft_save = fourier_cosine_transform (fr3, r, mshr, dr, k[0]);
		
		for (int ik = 0; ik < mshk-1; ik++)
		{
			double ft0 = ft_save;
			double dft0 = dft_save;
			
			double ft1 = fourier_sine_transform (fr2, r, mshr, dr, k[ik+1]);
			double dft1 = fourier_cosine_transform (fr3, r, mshr, dr, k[ik+1]);

			double c0, c1, c2, c3;
			c0 = ft0;
			c1 = dft0;
			c2 = 3.0*(ft1-ft0)/dk/dk-(dft1+2.0*dft0)/dk;
			c3 = (-2.0*(ft1-ft0)/dk+(dft1+dft0))/dk/dk;
			
			double k2 = k[ik]*k[ik];
			double k3 = k2*k[ik];
					
			polint_coef[ct] = c0-c1*k[ik]+c2*k2-c3*k3;
			polint_coef[ct+1] = c1-2.0*c2*k[ik]+3.0*c3*k2;
			polint_coef[ct+2] = c2-3.0*k[ik]*c3;
			polint_coef[ct+3] = c3;
			
//			double x = (k[ik]+k[ik+1])/2;
//			double x = k[ik];
//			double tmp = polint_coef[ct]+x*polint_coef[ct+1]+x*x*polint_coef[ct+2]+x*x*x*polint_coef[ct+3];
//			cout << "\ninterp = " << tmp << " ana = " << fourier_sine_transform (fr2, r, mshr, dr, x);
			//update
			ct += (polint_order+1);
			ft_save = ft1;
			dft_save = dft1;
		}
		
		//store coefficients for calculation
		double* coef = new double[n+1];
		double fac = dualfac (l) / dualfac (l-1);

		for (int j = 0; j < n; j++)
		{ 
			coef[j] = fac;

			//update
			int twoj = 2*j;
			fac *=  -static_cast<double>((l+twoj+2)*(l-twoj-1))/(twoj+3)/(twoj+2);

			//test
//			cout << "\ncoef[j] = " << coef[j] << endl;
		}
		coef[n] = pow (-1.0, n) * dualfac (2*l-1) / factorial (l);
		
		//start calc
		//special case k =0 ;
		fk[0] = 0.0;

		//k > 0
		for (int j = 0; j <= n; j++)
		{
			double Snm = 0.0;
			for (int ik = 1; ik < mshk; ik++)
			{
				Snm += pol_seg_int (polint_order, polint_coef, 2*j+1, k, ik);
				double k2j = pow(k[ik], 2*j+2);
				
				fk[ik] += coef[j] * Snm / k2j;
			}
		}
		
		//free
		delete [] coef;
	}

	delete [] fr2;
	delete [] fr3;
	delete [] polint_coef;
	
	timer::tick ("Mathzone_Add1","Sbt_new");
	return;
}

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

double Mathzone_Add1::pol_seg_int
(
	const int& polint_order,
	const double* coef,
	const int& n,
	const double* k,
	const int& ik
)
{
	double val = 0.0;
	double kmf = pow (k[ik], n+1);
	double kmb = pow (k[ik-1], n+1);
	
	int cstart = (polint_order+1)*(ik-1);
	for (int i = 0; i <= polint_order; i++)
	{
		val += coef[cstart+i]*(kmf-kmb)/(n+i+1);
		/*
		if (ik == 110) 
		{
			cout << "i = " << i << " coef = " << coef[cstart+i] << " df = " << kmf-kmb << endl;
		}
		*/	
		//update
		kmf *= k[ik];
		kmb *= k[ik-1];
	}

	/*
	if (ik == 110) 
	{
		cout << "val = " << val << endl;
		QUIT ();
	}
	*/
	return val;
}

double Mathzone_Add1::fourier_sine_transform
(
	const double* func,
	const double* r,
	const int& mshr,
	const double& dr,
	const double& k
)
{
	timer::tick ("Mathzone_Add1","Fsin");
	double val = 0.0;
	double* sinf = new double[mshr];
	for (int ir = 0; ir < mshr; ir++)
	{
		sinf[ir] = func[ir] * sin(k*r[ir]);
	}
	val = uni_simpson (sinf, mshr, dr);	
	delete[] sinf;
	timer::tick ("Mathzone_Add1","Fsin");
	return val;
}

double Mathzone_Add1::fourier_cosine_transform
(
	const double* func,
	const double* r,
	const int& mshr,
	const double& dr,
	const double& k
)
{
	timer::tick ("Mathzone_Add1","Fcos");
	double val = 0.0;
	double* cosf = new double[mshr];
	for (int ir = 0; ir < mshr; ir++)
	{
		cosf[ir] = func[ir] * cos(k*r[ir]);
	}

	val = uni_simpson (cosf, mshr, dr);	
	delete[] cosf;
	timer::tick ("Mathzone_Add1","Fcos");
	return val;
}

void Mathzone_Add1::test ()
{
	int polint_order =3;
	int dim = 2048;
	int ci = 1;
	int l = 0;
	double rmax = 20;
	double dr = rmax/dim;
	
	double* rad = new double[dim];
	double* func = new double[dim];
	double* fk = new double[dim];
	for (int ir = 0; ir < dim; ir++)
	{
		rad[ir] = ir * dr;
		func[ir] = pow(rad[ir], l) * exp(-ci*rad[ir]*rad[ir]);
		fk[ir] = 0.0;
	}
	
	Sbt_new (polint_order, l, rad, dr, dim, rad, dr, dim, func, 0, fk);
	
	for (int ik = 0; ik < dim; ik++)
	{
		double diff = fk[ik]- sqrt(PI/4/ci)/pow(2.0*ci, l+1)* std::pow(rad[ik], l) * exp(-rad[ik]*rad[ik]/4/ci);
		cout << rad[ik] << " " << fk[ik] << " " << sqrt(PI/4/ci)/pow(2.0*ci, l+1)*pow(rad[ik], l)*exp(-rad[ik]*rad[ik]/4/ci)
		<< " "	<< std::log(fabs(diff))/std::log(10.0) << endl;
	}

	delete[] rad;
	delete[] func;
	delete[] fk;
	return;
}

void Mathzone_Add1::test2 ()
{
	int polint_order =3;
	int N = 200;
	int ci = 1;
	int l = 0;
	double rmax = 20;
	double dr = rmax/(N-1);

	double dk = PI / rmax /2;
//	double kmax = PI / dr; 
//	double dk = dr;
	
	double* rad = new double[N];
	double* kad = new double[N];
	double* func = new double[N];
	double* fk = new double[N];
	double* fr = new double[N];
	for (int ir = 0; ir < N; ir++)
	{
		rad[ir] = ir * dr;
		kad[ir] = ir * dk;
		func[ir] = pow(rad[ir], l) * exp(-ci*rad[ir]*rad[ir]);
		fk[ir] = 0.0;
		fr[ir] = 0.0;
	}
	
	Sbt_new (polint_order, l, kad, dk, N, rad, dr, N, func, 0, fk);
	
/*	
	for (int ik = 0; ik < N; ik++)
	{
		double diff = fk[ik]- sqrt(PI/4/ci)/pow(2*ci, l+1)*pow(kad[ik], l)*exp(-kad[ik]*kad[ik]/4/ci);
		cout << kad[ik] << " " << fk[ik] << " " << sqrt(PI/4/ci)/pow(2*ci, l+1)*pow(kad[ik], l)*exp(-kad[ik]*kad[ik]/4/ci)
		<< " "	<< log(fabs(diff))/log(10) << endl;
	}
	QUIT ();
*/
	Sbt_new (polint_order, l, rad, dr, N, kad, dk, N, fk, 0, fr);

	for (int ir = 0; ir < N; ir++)
	{
		cout << ir*dr << " " << func[ir] << " " << fr[ir] *2.0 / PI << " " << std::log(fabs(fr[ir]*2.0/PI-func[ir]))/std::log(10.0) << endl;
	}
	

	delete[] rad;
	delete[] func;
	delete[] fk;
	delete[] fr;
	delete[] kad;
	return;
}

double Mathzone_Add1::Polynomial_Interpolation
(
 	const double* xa,
	const double* ya,
	const int& n,
	const double& x
)
{
	timer::tick("Mathzone_Add1","Polynomial_Interpolation");

	int i, m, ns;
	double 	den, dif, dift, ho, hp, w, rs, drs;
	
	//zero offset
	const double* Cxa = xa - 1;
	const double* Cya = ya - 1;
	
	double* cn = new double[n+1];
	double* dn = new double[n+1];

	ns = 1;
	dif = fabs(x - Cxa[1]);
	
	for(i = 1; i <= n; i++)
	{
		dift = fabs(x - Cxa[i]);
		if(dift < dif)
		{
			ns = i;
			dif = dift;
		}
		cn[i] = Cya[i];
		dn[i] = Cya[i];
	}

	rs = Cya[ns--];

	for(m = 1; m < n; m++)
	{
		for(i = 1; i <= n-m; i++)
		{
			ho = Cxa[i] - x;
			hp = Cxa[i+m] - x;
			w = cn[i+1] - dn[i];

			den = ho - hp;
			if(den == 0.0) 
			{
				cout << "Two Xs are equal" << endl;
				// WARNING_QUIT("Mathzone_Add1::Polynomial_Interpolation","Two Xs are equal");
				exit(0); // mohan update 2021-05-06
			}
			den = w / den;

			dn[i] = hp * den;
			cn[i] = ho * den;
		}
		if(2 * ns < n-m) drs = cn[ns+1];
		else drs = dn[ns--];

		rs += drs;
	}

	delete[] cn;
	delete[] dn;
	
	return rs;
	timer::tick("Mathzone_Add1","Polynomial_Interpolation");

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
	timer::tick("Mathzone_Add1","SplineD2");
	
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

	timer::tick("Mathzone_Add1","SplineD2");
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
	timer::tick("Mathzone_Add1","Cubic_Spline_Interpolation");

	#pragma omp parallel for schedule(static)
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
                
		//if(klo==0 && khi==1)                    // extrapolation
		//{ 
		//        klo = 1;
		//        khi = 2;
		//}

		const double h = rad[khi] - rad[klo];
		if(h == 0.0) 
		{
			cout << "Cubic_Spline_Interpolation, h == 0.0 so that cannot be divided" << endl;
			//WARNING_QUIT("Cubic_Spline_Interpolation","h == 0.0 so that cannot be divided");
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
	}

	timer::tick("Mathzone_Add1","Cubic_Spline_Interpolation");
}

// Interpolation for Numerical Orbitals
double Mathzone_Add1::RadialF
(
	const double* rad,
	const double* rad_f,
	const int& msh,
	const int& l,
	const double& R
)
{
  timer::tick("Mathzone_Add1","RadialF");
  
  int mp_min, mp_max, m;
  double h1, h2, h3, f1, f2, f3, f4;
  double g1, g2, x1, x2, y1, y2, f;
  double c, result;

  mp_min = 0;
  mp_max = msh - 1;

  //assume psir behaves like r**l
  if (R < rad[0])
  {
    if (l == 0) 
	{
		f = rad_f[0];
	}
	else 
	{
		c = rad_f[0] / pow(rad[0], l);
		f = pow(R, l) * c;
	}
  }
  else if (rad[mp_max] < R) 
  {
	  f = 0.0;
  }
  else
  {
	  do
	  {
		  m = (mp_min + mp_max)/2;
		  if (rad[m] < R) mp_min = m;
		  else mp_max = m;
	  }
	  while((mp_max-mp_min)!=1);
	  m = mp_max;

	  if (m < 2) 
	  {
		  m = 2;
	  }
	  else if (msh <= m) 
	  {
		  m = msh - 2;
	  }

	  /****************************************************
		Spline like interpolation
	   ****************************************************/

	  if (m == 1)
	  {
		  h2 = rad[m] - rad[m-1];
		  h3 = rad[m+1] - rad[m];

		  f2 = rad_f[m-1];
		  f3 = rad_f[m];
		  f4 = rad_f[m+1];

		  h1 = -(h2+h3);
		  f1 = f4;
	  }
	  else if (m == (msh-1))
	  {
		  h1 = rad[m-1] - rad[m-2];
		  h2 = rad[m] - rad[m-1];

		  f1 = rad_f[m-2];
		  f2 = rad_f[m-1];
		  f3 = rad_f[m];

		  h3 = -(h1+h2);
		  f4 = f1;
	  }
	  else
	  {
		  h1 = rad[m-1] - rad[m-2];
		  h2 = rad[m]   - rad[m-1];
		  h3 = rad[m+1] - rad[m];

		  f1 = rad_f[m-2];
		  f2 = rad_f[m-1];
		  f3 = rad_f[m];
		  f4 = rad_f[m+1];
	  }

	  //Calculate the value at R

	  g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
	  g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

	  x1 = R - rad[m-1];
	  x2 = R - rad[m];
	  y1 = x1/h2;
	  y2 = x2/h2;

	  f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
		  + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }

  result = f;

  timer::tick("Mathzone_Add1","RadialF");
  return result;
}

// Interpolation for Numerical Orbitals
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
		throw runtime_error("newr should >= 0. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));

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
	timer::tick("Mathzone_Add1", "Uni_Deriv_Phi");
	int FFT_NR = 2*mesh-1;  // FFT_NR = 16019
	// cout << "\n mesh=" << mesh << ", radf[8010]=" << radf[8010] <<  ", radf[8009]=" << radf[8009] ;
	// mesh=8010, radf[8010]=4.396478951532926e-01, radf[8009]=0.000000000000000e+00

	fftw_complex fft_phir[FFT_NR], fft_phik[FFT_NR];
	fftw_complex fft_ndphik[FFT_NR], fft_ndphir[FFT_NR];
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
	//cout << "\n Call FFTW3 ";
	p1 = fftw_plan_dft_1d(FFT_NR, fft_phir, fft_phik, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p1);
	//fftw_destroy_plan(p1);
#elif defined __FFTW2
	//cout << "\n Call FFTW2 ";
	p1 = fftw_create_plan(FFT_NR, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_one(p1, fft_phir, fft_phik);
	//fftw_destroy_plan(p1);
#endif

	
	double dk_uniform = TWO_PI / FFT_NR / dr;
	
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
	
	timer::tick("Mathzone_Add1", "Uni_Deriv_Phi");
}
