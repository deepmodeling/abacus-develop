#include "integral.h"
#include "polint.h"
#include "constants.h"
#include <cassert>
#include <cmath>
using namespace std;

int Integral::n_root = 512;
bool Integral::calc_wx = false;
double* Integral::gauleg_w;
double* Integral::gauleg_x;

double Integral::Gauss_Legendre
(
	const double &a,
	const double &b,
	const double* F,
	const double* Rad,
	const int& Msh
)
{
	assert(Msh > 0);
	
	int ir;
	double sab, dab, xgl, ygl, sum;
	const double tiny = 1e-9;
	
	if(!calc_wx) 
	{
		Integral::gauleg();
		calc_wx = true;
	}

	sab = a + b;
	dab = b - a;
	//test
//	for(ir = 0; ir < n_root; ir++)
//	{
//		cout << gauleg_x[ir] << " " << gauleg_w[ir] << endl;
//	}
	if( fabs(sab) < tiny ) return 0.0;

	//initialization
	sum = 0.0;
	for(ir = 0; ir < n_root; ir++)
	{
		xgl = (sab + dab * gauleg_x[ir]) / 2;
		ygl = Polint::RadialF(Rad, F, Msh, 0, xgl);
		ygl = Polint::Lagrange3(Rad, F, Msh, xgl);
//		cout << "\nxgl = " << xgl << " ygl = " << ygl << endl;

		sum += gauleg_w[ir] * ygl;
	}

	return sum * dab / 2.0;
}

void Integral::gauleg()
{
	int m, j, i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	const double x1 = -1.0;
	const double x2 = 1.0;
	const double EPS = 1.0e-14;

	gauleg_x = new double[n_root+1];
	gauleg_w = new double[n_root+1];
	m=(n_root + 1) / 2;
	xm=0.50*(x2+x1);
	xl=0.50*(x2-x1);

	for(i = 1; i <= m; i++)
	{
		z = cos( PI * (i - 0.250) / (n_root + 0.5) );

		do
		{
			p1 = 1.0;
			p2 = 0.0;

			for( j = 1; j <= n_root; j++)
			{
				p3 = p2;
				p2 = p1;
				p1=( (2.0 * static_cast<double>(j) - 1.0) * z * p2 -
						( static_cast<double>(j) - 1.0 ) *p3  ) / static_cast<double>(j);
			}
			pp = static_cast<double>( n_root ) * ( z * p1 - p2 ) / ( z * z - 1.0 );
			z1 = z;
			z = z1 - p1 / pp;
		} while( fabs(z - z1) > EPS );

		gauleg_x[i] = xm - xl * z;
		gauleg_x[n_root + 1 - i] = xm + xl * z;
		gauleg_w[i] = 2.0 * xl / ( (1.0 - z*z) * pp *pp );
		gauleg_w[n_root + 1 - i] = gauleg_w[i];

	}

	for(i = 1; i <= n_root; i++)
	{
		gauleg_x[i-1] = gauleg_x[i];
		gauleg_w[i-1] = gauleg_w[i];
	}
	return;
}

