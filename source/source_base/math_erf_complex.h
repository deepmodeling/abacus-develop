#ifndef MATH_ERF_COMPLEX_H
#define MATH_ERF_COMPLEX_H

#include <complex>
#include <vector>

namespace ModuleBase
{

class ErrorFunc
{
  public:
    ErrorFunc();
    ~ErrorFunc();

    /**
     * @brief A class of the error function of complex arguments based on Faddeeva algorithm.
     *       More information please refer to http://ab-initio.mit.edu/Faddeeva
     *
     * @author jiyy on 2025-05-18
     */

    // compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
    static std::complex<double> scaled_w(std::complex<double> z, double relerr);
    static double scaled_w_im(double x); // special-case code for Im[w(x)] of real x

    // compute erfcx(z) = exp(z^2) erfc(z)
    static std::complex<double> erfcx(std::complex<double> z, double relerr = 0);
    static double erfcx(double x); // special case for real x

    // compute erf(z), the error function of complex arguments
    static std::complex<double> erf(std::complex<double> z, double relerr = 0);

    // compute erfi(z) = -i erf(iz), the imaginary error function
    static std::complex<double> erfi(std::complex<double> z, double relerr = 0);
    static double erfi(double x); // special case for real x

    // compute erfc(z) = 1 - erf(z), the complementary error function
    static std::complex<double> erfc(std::complex<double> z, double relerr = 0);

  private:
    static double w_im_y100(double y100, double x);
    static inline double sinc(double x, double sinx)
    {
        return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667) * x * x : sinx / x;
    }
    static inline double copysign(double x, double y)
    {
        return x < 0 != y < 0 ? -x : x;
    }
    static const std::vector<double> expa2n2;
};
} // namespace ModuleBase

#endif // MATH_ERF_H
