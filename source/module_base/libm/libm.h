//==========================================================
// AUTHOR : alcanderian@gmail.com
// DATE : 2023-01-06
//==========================================================

#ifndef MODULE_BASE_ABAUCS_LIBM_H
#define MODULE_BASE_ABAUCS_LIBM_H

#include <complex>
#include <cmath>

namespace ModuleBase
{
namespace libm
{

#ifdef USE_ABACUS_LIBM

double __exp (double x);
double __cos (double x);
double __sin (double x);
void __sincos (double x, double *sinx, double *cosx);
std::complex<double> __cexp (const std::complex<double> &x);

#else

inline double __exp (double x) { return std::exp(x); };
inline double __cos (double x) { return std::cos(x); };
inline double __sin (double x) { return std::sin(x); };
inline void __sincos (double x, double *sinx, double *cosx) { sincos(x, sinx, cosx); };
inline std::complex<double> __cexp (const std::complex<double> &x) { return std::exp(x); }

#endif

inline float __expf (float x) { return std::exp(x); };
inline float __cosf (float x) { return std::cos(x); };
inline float __sinf (float x) { return std::sin(x); };
inline void __sincosf (float x, float *sinx, float *cosx) { sincosf(x, sinx, cosx); };
inline std::complex<float> __cexpf (const std::complex<float> &x) { return std::exp(x); }

template<typename Tp> Tp exp(Tp x);
template<typename Tp> Tp cos(Tp x);
template<typename Tp> Tp sin(Tp x);
template<typename Tp> void sincos(Tp x, Tp *s, Tp *c);
template<typename Tp> std::complex<Tp> exp(const std::complex<Tp> &x);

template<> inline double exp(double x) { return __exp(x); }
template<> inline double cos(double x) { return __cos(x); }
template<> inline double sin(double x) { return __sin(x); }
template<> inline void sincos(double x, double *s, double *c) { __sincos(x, s, c); }
template<> inline std::complex<double> exp(const std::complex<double> &x) { return __cexp(x); }

template<> inline float exp(float x) { return __expf(x); }
template<> inline float cos(float x) { return __cosf(x); }
template<> inline float sin(float x) { return __sinf(x); }
template<> inline void sincos(float x, float *s, float *c) { __sincosf(x, s, c); }
template<> inline std::complex<float> exp(const std::complex<float> &x) { return __cexpf(x); }

};
};

#endif
