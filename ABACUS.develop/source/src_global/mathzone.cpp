//#include "../src_pw/global.h"
#include "mathzone.h"
#include "../module_base/constants.h"
#include "../module_base/global_variable.h"
#include "../module_base/global_function.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"

Mathzone::Mathzone()
{
}

Mathzone::~Mathzone()
{
}

void Mathzone::norm_pw(complex<double> *u,const int n)
{
    if (n < 1) 
	{
		WARNING_QUIT("Mathzone::norm","n < 1");
	}

    complex<double> rrc = ZERO;//complex < double> (0.0,0.0);

    for (int i = 0;i < n;i++) 
	{
		rrc += conj(u[i]) * u[i];
	}

    Parallel_Reduce::reduce_complex_double_pool( rrc );

    double rr = rrc.real();

    if (rr <= 1.e-20) 
	{
		WARNING_QUIT("Mathzone::norm","rr <= 1.e-20");
	}

    rr = 1.0/sqrt(rr);

    if (rr==1.0) 
	{
		return;
	}
    else
    {
        for (int i = 0;i < n;i++) 
		{
			u[i] *= rr;
		}
    }
	return;
}


void Mathzone::To_Polar_Coordinate
(
    const double &x_cartesian,
    const double &y_cartesian,
    const double &z_cartesian,
    double &r,
    double &theta,
    double &phi
)
{
//	OUT("x",x_cartesian);
//	OUT("y",y_cartesian);
//	OUT("z",z_cartesian);
    r = sqrt(x_cartesian*x_cartesian
             + y_cartesian*y_cartesian
             + z_cartesian*z_cartesian);
//----------------------------------------------------------
// Calculate theta
//----------------------------------------------------------
    double small = 1.0e-9;
    if (r < small)
    {
        theta = PI/2.0;
    }
    else
    {
        theta = acos(z_cartesian / r);
    }
//----------------------------------------------------------
// Calculate phi
//----------------------------------------------------------
    if (x_cartesian > small && y_cartesian > small)
    {
        phi = atan( y_cartesian / x_cartesian );
    }
    else if (  x_cartesian < small  )
    {
        phi = atan( y_cartesian / x_cartesian ) + PI;
    }
    else if ( x_cartesian > small && y_cartesian < -small)
    {
        phi = atan( y_cartesian/x_cartesian ) + TWO_PI;
    }
    else
    {
        phi = PI/2.0 * ( (y_cartesian >= 0.0) ? 1.0 : -1.0);
    }
//----------------------------------------------------------
// degress = radians / PI * 180
//----------------------------------------------------------
    theta = theta/PI*180;
    phi = phi/PI*180;

//	OUT("r",r);
//	OUT("theta(degree)",theta);
//	OUT("phi(degree)",phi);
    return;
}

