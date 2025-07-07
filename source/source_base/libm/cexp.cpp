//==========================================================
// AUTHOR : alcanderian@gmail.com
// DATE : 2023-01-06
//==========================================================

#include <complex>
#include <math.h>
#include <stdint.h>
#include <float.h>
#include <limits>

namespace ModuleBase
{
namespace libm
{

template<typename FLOAT> void __msincos (FLOAT x, FLOAT *s, FLOAT *c);
template<typename FLOAT> FLOAT __mexp (FLOAT x);
template<typename FLOAT> FLOAT __mhugeval ();

void __sincos (double x, double *sinx, double *cosx);
double __exp (double x);

template<> inline double __mexp (double x) { return __exp (x); }
template<> inline void __msincos (double x, double *s, double *c) { __sincos (x, s, c); }
template<> inline double __mhugeval () { return HUGE_VAL; }

template<typename FLOAT>
inline std::complex<FLOAT>
__cexp_impl (const std::complex<FLOAT> &x)
{
  std::complex<FLOAT> retval;
  int rcls = fpclassify (x.real());
  int icls = fpclassify (x.imag());

  if (rcls >= FP_ZERO)
    {
      /* Real part is finite.  */
      if (icls >= FP_ZERO)
	{
      
	  /* Imaginary part is finite.  */
	  const int t = (int) ((std::numeric_limits<FLOAT>::max_exponent - 1) * (FLOAT)(M_LN2));
	  FLOAT sinix, cosix;

	  if ((fabs (x.imag()) > std::numeric_limits<FLOAT>::min()))
	    {
	      __msincos (x.imag(), &sinix, &cosix);
	    }
	  else
	    {
	      sinix = x.imag();
	      cosix = 1;
	    }

      std::complex<FLOAT> __x = x;
	  if (__x.real() > t)
	    {
	      FLOAT exp_t = __mexp ((FLOAT)t);
	      __x -= std::complex<FLOAT>(t, 0);
	      sinix *= exp_t;
	      cosix *= exp_t;
	      if (__x.real() > t)
		{
		  __x -= std::complex<FLOAT>(t, 0);
		  sinix *= exp_t;
		  cosix *= exp_t;
		}
	    }
	  if (__x.real() > t)
	    {
	      /* Overflow (original real part of x > 3t).  */
	      retval.real(std::numeric_limits<FLOAT>::max() * cosix);
	      retval.imag(std::numeric_limits<FLOAT>::max() * sinix);
	    }
	  else
	    {
	      FLOAT exp_val = __mexp (__x.real());
	      retval.real(exp_val * cosix);
	      retval.imag(exp_val * sinix);
	    }
	}
      else
	{
	  /* If the imaginary part is +-inf or NaN and the real part
	     is not +-inf the result is NaN + iNaN.  */
	  retval.real(std::numeric_limits<FLOAT>::quiet_NaN());
	  retval.imag(std::numeric_limits<FLOAT>::quiet_NaN());
	}
    }
  else if (rcls == FP_INFINITE)
    {
      /* Real part is infinite.  */
      if (icls >= FP_ZERO)
	{
	  /* Imaginary part is finite.  */
	  FLOAT value = signbit (x.real()) ? 0 : __mhugeval<FLOAT>();

	  if (icls == FP_ZERO)
	    {
	      /* Imaginary part is 0.0.  */
	      retval.real(value);
	      retval.imag(x.imag());
	    }
	  else
	    {
	      FLOAT sinix, cosix;

	      if ((fabs (x.imag()) > std::numeric_limits<FLOAT>::min()))
		{
		  __msincos (x.imag(), &sinix, &cosix);
		}
	      else
		{
		  sinix = x.imag();
		  cosix = 1;
		}

	      retval.real(copysign (value, cosix));
	      retval.imag(copysign (value, sinix));
	    }
	}
      else if (signbit (x.real()) == 0)
	{
	  retval.real(__mhugeval<FLOAT>());
	  retval.imag(x.imag() - x.imag());
	}
      else
	{
	  retval.real(0);
	  retval.imag(copysign (0, x.imag()));
	}
    }
  else
    {
      /* If the real part is NaN the result is NaN + iNaN unless the
	 imaginary part is zero.  */
      retval.real(std::numeric_limits<FLOAT>::quiet_NaN());
      if (icls == FP_ZERO)
	retval.imag(x.imag());
      else
	{
	  retval.imag(std::numeric_limits<FLOAT>::quiet_NaN());
	}
    }

  return retval;
}

std::complex<double>
__cexp (const std::complex<double> &x) { return __cexp_impl(x); }

};
};
