//==========================================================
// AUTHOR : Peize Lin
// DATE : 2017-04-24
//==========================================================
#include "sph_bessel_recursive.h"

#include <cmath>
#include <stdexcept>

#include "constants.h"

namespace ModuleBase
{

std::vector<Sph_Bessel_Recursive::D1> Sph_Bessel_Recursive_Pool::D1::sb_pool;

void Sph_Bessel_Recursive::D1::set_dx( const double dx_in )
{
	if(finish_set_dx && dx_in!=dx)	
		throw std::runtime_error("Sph_Bessel_Recursive::set_dx, dx can only set once");
	else
	{
		dx = dx_in;
		finish_set_dx = true;
	}
}

const std::vector<std::vector<double>> & Sph_Bessel_Recursive::D1::cal_jlx( const int lmax, const size_t ix_size )
{
	if(lmax<0)
		throw std::invalid_argument("Sph_Bessel_Recursive::jlx l<0");
	cal_jlx_0( lmax+1 );
	cal_jlx_smallx( lmax+1, ix_size );
	cal_jlx_recursive( lmax+1, ix_size );
	return jlx;
}

void Sph_Bessel_Recursive::D1::cal_jlx_0( const int l_size )
{
	if(jlx.size() < static_cast<size_t>(l_size))
		jlx.resize(l_size);
	
	for( int l=0; l!=l_size; ++l )
	{
		if(jlx[l].size()<1)
		{
			jlx[l].resize(1);
			if(0==l)
				jlx[l][0] = 1.0;
			else
				jlx[l][0] = 0.0;
		}
	}
}

void Sph_Bessel_Recursive::D1::cal_jlx_smallx( const int l_size, const size_t ix_size )
{
	if(jlx.size() < static_cast<size_t>(l_size))
		jlx.resize(l_size);
	
	for( int l=0; l!=l_size; ++l )
	{
		if(jlx[l].size()<ix_size)
		{
			const double coeff = sqrt(ModuleBase::PI)/tgamma(l+1.5)/pow(2,l+1);
			const double smallx_range = pow( this->threshold/coeff*(l+1.5)*4, 1.0/(l+2) );
			
			const size_t ix_size_begin = static_cast<int>(jlx[l].size());
			const size_t ix_size_end = std::min( ix_size, static_cast<size_t>(smallx_range/dx) );
			if(jlx[l].size()<ix_size_end)
			{
				jlx[l].resize(ix_size_end);
				for( size_t ix=ix_size_begin; ix<ix_size_end; ++ix )
				{
					const double x1 = ix * dx;
					jlx[l][ix] = coeff*pow(x1,l);
				}
			}
		}
	}	
}

void Sph_Bessel_Recursive::D1::cal_jlx_recursive( const int l_size, const size_t ix_size )
{
	if(jlx.size() < static_cast<size_t>(l_size))
		jlx.resize(l_size);
	
	for( int l=0; l!=l_size; ++l )
	{
		if(jlx[l].size()<ix_size)
		{
			const size_t ix_size_begin = static_cast<int>(jlx[l].size());
			jlx[l].resize(ix_size);
			switch(l)
			{
				case 0:
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						jlx[l][ix] = sin(x1)/x1;
					}
					break;
				case 1:
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						jlx[l][ix] = sin(x1)/x2-cos(x1)/x1;
					}
					break;
				default:
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						jlx[l][ix] = (2*l-1)/x1*jlx[l-1][ix] - jlx[l-2][ix];
					}
					break;
			}
		}
	}
}

}

/*
void Sph_Bessel_Recursive::cal_jlx_preset(const int l_size, const size_t ix_size)
{	
std::cout<<l_size<<"\t"<<ix_size<<std::endl;
	if(jlx.size()<l_size)
		jlx.resize(l_size);
	
	for( int l=0; l!=l_size; ++l )
	{
//timeval t_start; gettimeofday( &t_start, NULL);
		if(jlx[l].size()<ix_size)
		{
			const size_t ix_size_begin = static_cast<int>(jlx[l].size());
			jlx[l].resize(ix_size);
			switch(l)
			{
				case 0:
				{
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						jlx[l][ix] = sin(x1)*( 1.0/x1 );
					}
					break;
				}
				case 1:
				{
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						jlx[l][ix] = sin(x1)*( 1.0/x2 ) 
						           + cos(x1)*( - 1.0/x1 );
					}
					break;
				}
				case 2:
				{
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;						
						jlx[l][ix] = sin(x1)*(-1.0/x1+3.0/x3) 
						           + cos(x1)*(-3.0/x2);
					}
					break;
				}
				case 3:
				{
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						jlx[l][ix] = sin(x1)*( - 6.0/x2 + 15.0/x4 )
						           + cos(x1)*( 1.0/x1 - 15.0/x3 );
					}
					break;
				}
				case 4:
				{
					for( size_t ix=ix_size_begin; ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						jlx[l][ix] = sin(x1)*( 1.0/x1 - 45.0/x3 + 105.0/x5 )
						           + cos(x1)*( 10.0/x2 - 105.0/x4 );
					}
					break;
				}
				case 5:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.14/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						jlx[l][ix] = sin(x1)*( 15.0/x2 - 420.0/x4 + 945.0/x6 )
						           + cos(x1)*( - 1.0/x1 + 105.0/x3 - 945.0/x5 );
					}
					break;
				}
				case 6:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						jlx[l][ix] = sin(x1)*( - 1.0/x1 + 210.0/x3 - 4725.0/x5 + 10395.0/x7 )
						           + cos(x1)*( - 21.0/x2 + 1260.0/x4 - 10395.0/x6 );
					}
					break;
				}
				case 7:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						jlx[l][ix] = sin(x1)*( 135135.0/x8 - 28.0/x2 + 3150.0/x4 - 62370.0/x6 )
						           + cos(x1)*( 1.0/x1 - 378.0/x3 + 17325.0/x5 - 135135.0/x7 );
					}
					break;
				}
				case 8:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						const double x9 = x1 * x8;
						jlx[l][ix] = sin(x1)*( 2027025.0/x9 - 630.0/x3 + 1.0/x1 + 51975.0/x5 - 945945.0/x7 )
						           + cos(x1)*( - 2027025.0/x8 + 36.0/x2 - 6930.0/x4 + 270270.0/x6 );
					}
					break;
				}
				case 9:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						const double x9 = x1 * x8;
						const double x10 = x1 * x9;
						jlx[l][ix] = sin(x1)*( - 16216200.0/x8 + 34459425.0/x10 - 13860.0/x4 + 45.0/x2 + 945945.0/x6 )
						           + cos(x1)*( - 34459425.0/x9 + 990.0/x3 - 1.0/x1 - 135135.0/x5 + 4729725.0/x7 );
					}
					break;
				}
				case 10:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						const double x9 = x1 * x8;
						const double x10 = x1 * x9;
						const double x11 = x1 * x10;
						jlx[l][ix] = sin(x1)*( - 1.0/x1 + 1485.0/x3 - 315315.0/x5 + 18918900.0/x7 - 310134825.0/x9 + 654729075.0/x11 )
						           + cos(x1)*( 91891800.0/x8 - 654729075.0/x10 + 25740.0/x4 - 55.0/x2 - 2837835.0/x6 );
					}
					break;
				}
				case 11:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						const double x9 = x1 * x8;
						const double x10 = x1 * x9;
						const double x11 = x1 * x10;
						const double x12 = x1 * x11;
						jlx[l][ix] = sin(x1)*( - 66.0/x2 + 45045.0/x4 - 7567560.0/x6 + 413513100.0/x8 - 6547290750.0/x10 + 13749310575.0/x12 )
						           + cos(x1)*( 1.0/x1 - 2145.0/x3 + 675675.0/x5 - 64324260.0/x7 + 1964187225.0/x9 - 13749310575.0/x11 );
					}
					break;
				}
				case 12:
				{
					const size_t ix_size_piecewise = min( ix_size, static_cast<size_t>(0.29/dx)+1 );
					for( size_t ix=ix_size_begin; ix<ix_size_piecewise; ++ix )
						jlx[l][ix] = 0;			//mohan add 2007-10-15
					for( size_t ix=max(ix_size_begin,ix_size_piecewise); ix<ix_size; ++ix )
					{
						const double x1 = ix * dx;
						const double x2 = x1 * x1;
						const double x3 = x1 * x2;
						const double x4 = x1 * x3;
						const double x5 = x1 * x4;
						const double x6 = x1 * x5;
						const double x7 = x1 * x6;
						const double x8 = x1 * x7;
						const double x9 = x1 * x8;
						const double x10 = x1 * x9;
						const double x11 = x1 * x10;
						const double x12 = x1 * x11;
						const double x13 = x1 * x12;
//						jlx[l][ix] = sin(x1)*( 1.0/x1 - 3003.0/x3 + 1351350.0/x5 - 192972780.0/x7 + 9820936125.0/x9 - 151242416325.0/x11 + 316234143225.0/x13 )
//						           + cos(x1)*( 78.0/x2 - 75075.0/x4 + 18378360.0/x6 - 1571349780.0/x8 + 45831035250.0/x10 - 316234143225.0/x12 );
						jlx[l][ix] = sin(x1)*( (1.0 +(-3003.0 +(1351350.0/x5 +(-192972780.0 +(9820936125.0 +(-151242416325.0 +316234143225.0/x2)/x2)/x2)/x2)/x2)/x2)/x1 )
						           + cos(x1)*( (78.0 +(-75075.0 +(18378360.0 +(-1571349780.0 +(45831035250.0 -316234143225.0/x2)/x2)/x2)/x2)/x2)/x2 );
					}
					break;
				}
				default:
					throw std::invalid_argument("Sph_Bessel_Recursive::jlx l<0 or l>6");
			}
		}
//timeval t_end; gettimeofday( &t_end, NULL);
//std::cout<<l<<"\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
	}
}
*/	
