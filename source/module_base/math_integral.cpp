#include "math_integral.h"
#include <stddef.h> // use size_t
#include <cassert>
#include <algorithm>
#include <functional>
#include <cmath>

namespace ModuleBase
{

Integral::Integral(){}

Integral::~Integral(){}


// Peize Lin accelerate 2017-10-02
/*
void Integral::Simpson_Integral
(
    const int mesh,
    const double *func,
    const double *rab,
    double &asum
)
{
    //     simpson's rule integration. On input:
    //     mesh = mhe number of grid points (should be odd)
    //     func(i)= function to be integrated
    //     rab(i) = r(i) * dr(i)/di * di
    //     For the logarithmic grid not including r=0 :
    //     r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    //     For the logarithmic grid including r=0 :
    //     r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    //     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    //     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    
    //  simpson's rule integrator for function stored on the
    //  radial logarithmic mesh
    //	routine assumes that mesh is an odd number so run check
    if (mesh % 2 == 0)
    {
        std::cout << "\n error in subroutine simpson ";
        std::cout << "\n routine assumes mesh is odd but mesh = "
             << mesh << std::endl;
        return;
    }

    asum = 0.00;
    const double r12 = 1.00 / 12.00;
    double f3 = func [0] * rab [0] * r12;
    for (int i = 1;i < mesh;i += 2)
    {
        const double f1 = f3;
        const double f2 = func [i] * rab [i] * r12;
        f3 = func [i + 1] * rab [i + 1] * r12;
        asum += 4.00 * f1 + 16.00 * f2 + 4.00 * f3;
    }
    return;
}// end subroutine simpson
*/


// Peize Lin accelerate 2017-10-02
void Integral::Simpson_Integral
(
    const int mesh,
    const double * const func,
    const double * const rab,
    double &asum
)
{
    /*     simpson's rule integration. On input:
    !      mesh = mhe number of grid points (should be odd)
    !      func(i)= function to be integrated
    !      rab(i) = r(i) * dr(i)/di * di
    !      For the logarithmic grid not including r=0 :
    !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    !      For the logarithmic grid including r=0 :
    !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    */
    //  simpson's rule integrator for function stored on the
    //  radial logarithmic mesh
    //	routine assumes that mesh is an odd number so run check
    assert(mesh&1);

    asum = 0.00;
	const size_t end = mesh-2;
    for( size_t i=1; i!=end; i+=2 )
    {
		const double f1 = func[i]*rab[i];
		asum += f1 + f1 + func[i+1]*rab[i+1];
    }
	const double f1 = func[mesh-2]*rab[mesh-2];
	asum += f1+f1;
	asum += asum;
	asum += func[0]*rab[0] + func[mesh-1]*rab[mesh-1];
	asum /= 3.0;
    return;
}// end subroutine simpson


// Peize Lin accelerate 2017-10-02
void Integral::Simpson_Integral
(
    const int mesh,
    const double * const func,
    const double dr,
    double &asum
)
{
    /*     simpson's rule integration. On input:
    !      mesh = mhe number of grid points (should be odd)
    !      func(i)= function to be integrated
    !      rab(i) = r(i) * dr(i)/di * di
    !      For the logarithmic grid not including r=0 :
    !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
    !      For the logarithmic grid including r=0 :
    !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
    !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
    !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
    */
    //  simpson's rule integrator for function stored on the
    //  radial logarithmic mesh
    //	routine assumes that mesh is an odd number so run check
    assert(mesh&1);

    asum = 0.00;
	const size_t end = mesh-2;
    for(size_t i=1; i!=end; i+=2 )
    {
		const double f1 = func[i];
		asum += f1 + f1 + func[i+1];
    }
	const double f1 = func[mesh-2];
	asum += f1+f1;
	asum += asum;
	asum += func[0] + func[mesh-1];
	asum *= dr/3.0;
    return;
}// end subroutine simpson


// Peize Lin add 2016-02-14
void Integral::Simpson_Integral_0toall
(
    const int mesh,
    const double * const func,
    const double * const rab,
    double * const asum
)
{
    // asum(r) = \int_{r'=0}^{r} dr' f(r') 

    const double r2=1.00/2.00, r3=1.00/3.00;
    asum[0] = 0.00;
    double f3 = func [0] * rab [0];
    for( int i=1; i<mesh; i+=2)
    {
        const double f1 = f3;
        const double f2 = func[i] * rab[i] ;
        f3 = func[i+1] * rab[i+1] ;
        asum[i] = asum[i-1] + r2*( f1 + f2);
        if(i+1<mesh)
        {
            asum[i+1] = asum[i-1] + r3*( f1 + 4.00*f2 + f3 );
        }
    }
    return;
}


// Peize Lin add 2016-02-14
// faster but still have bug
/*void Integral::Simpson_Integral_alltoinf
(
    const int mesh,
    const double *func,
    const double *rab,
    double *asum
)
{
    // asum(r) = \int_{r'=r}^{+\infty} dr' f(r') 
    //         = \inf_{r'=r}^{mesh} dr' f(r')

    const double r2=1.00/2.00, r3=1.00/3.00;
    asum[mesh-1] = 0.00;
    const int odd_mesh = (mesh-1)^~1;
    double f1 = func[odd_mesh] * rab[odd_mesh];
    for( size_t i=(mesh-3)|1; i>0; i-=2)
    {
        const double f3 = f1;   
        if( i+3==mesh )
        {
            const double f4 = func[mesh-1] * rab[mesh-1];
            asum[mesh-2] = r2*(f3 + f4);
        }
        const double f2 = func[i] * rab[i] ;
        f1 = func[i-1] * rab[i-1] ;
        asum[i-1] = asum[i+1] + r3*( f1 + 4.00*f2 + f3 );
        asum[i] = asum[i-1] - r2*( f1 + f2);
    }
    return;
}*/


// Peize Lin add 2016-06-11
// a little lower
void Integral::Simpson_Integral_alltoinf
(
    const int mesh,
    const double * const func,
    const double * const rab,
    double * const asum
)
{
    Integral::Simpson_Integral_0toall( mesh, func, rab, asum );

    const double asum_all = asum[mesh-1];
    for (int i = 0;i < mesh; ++i)
	{
        asum[i] = asum_all - asum[i];
	}
	return;
}

double Integral::simpson(const int n, const double* const f, const double dx) 
{
    assert(n >= 2);
    if (n == 4)
    { // Simpson's 3/8 rule
        return 3.0 * dx / 8 * (f[0] + 3.0 * f[1] + 3.0 * f[2] + f[3]);
    }

    if (n == 2)
    {
        return 0.5 * dx * (f[0] + f[1]);
    }

    if (n % 2 == 1)
    { // composite Simpson's 1/3 rule
        double sum = 0.0;
        for (int i = 1; i != n-2; i += 2) 
        {
            sum += 2.0 * f[i] + f[i+1];
        }
        sum += 2.0 * f[n-2];
        sum *= 2.0;
        sum += f[0] + f[n-1];
        return sum * dx / 3.0;
    }
    else
    { // composite Simpson's 1/3 rule for the first n-4 intervals plus Simpson's 3/8 rule for the last 3 intervals
        return simpson(n-3, f, dx) + simpson(4, &f[n-4], dx);
    }
}

double Integral::simpson(const int n, const double* const f, const double* const h)
{   
    // Simpson's rule for irregularly-spaced grid
    // The treatment for even number of grid points is the same as that of the regularly-spaced grid case.
    assert( n >= 2 );
    assert( std::all_of(h, h+(n-1), [](double x){return x > 0.0;}) );

    if (n == 4)
    {
        double w = h[0] + h[1] + h[2];
        return w / 12.0 * ( 2.0 + ((h[1]+h[2])/h[0]-1.0) * (h[2]/(h[0]+h[1])-1.0) ) * f[0]
            + std::pow(w,3) / 12.0 * (h[0]+h[1]-h[2]) / (h[0]*h[1]*(h[1]+h[2])) * f[1]
            + std::pow(w,3) / 12.0 * (h[2]+h[1]-h[0]) / (h[2]*h[1]*(h[1]+h[0])) * f[2]
            + w / 12.0 * ( 2.0 + ((h[1]+h[0])/h[2]-1.0) * (h[0]/(h[2]+h[1])-1.0) ) * f[3];
    }

    if (n == 2)
    {
        return 0.5 * h[0] * (f[0] + f[1]);
    }

    if (n % 2 == 1)
    {
        double sum = 0.0;
        for (int i = 0; i < n/2; ++i)
        {
            double hrp = h[2*i+1] / h[2*i];
            double hrm = h[2*i] / h[2*i+1];
            sum += (h[2*i+1] + h[2*i]) / 6.0 * ( (2.0-hrp)*f[2*i] + (2.0+hrp+hrm)*f[2*i+1] + (2.0-hrm) * f[2*i+2]);
        }
        return sum;
    }
    else
    {
        return simpson(n-3, f, h) + simpson(4, &f[n-4], &h[n-4]);
    }
}

}
