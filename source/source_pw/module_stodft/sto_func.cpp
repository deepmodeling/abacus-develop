#include "sto_func.h"

#include "source_estate/occupy.h"

#include <cmath>

#define TWOPI 6.283185307179586

namespace
{
// Evaluate the Fermi-Dirac entropy density  f*ln(f) + (1-f)*ln(1-f)
// from  ne_mu = (E - mu)/kT  without the cancellation that the literal
// expression suffers when f is tiny:
//     (1-f)*ln(1-f)  ~  -f,
// but in double precision 1-f rounds to 1 as soon as f < ~1e-16, so the
// literal expression drops this -f term and leaves an error of order f
// (relative error ~1/|ln f|, i.e. percent level) in the deep tail.
// The resulting non-smooth error seeds spurious high-order Chebyshev
// coefficients (~1e-12) which are then amplified by T_n(x) for any
// eigenvalue outside the expansion interval.  Using s = min(f, 1-f)
// computed directly (no subtraction) and log1p(-s) keeps full relative
// accuracy for all temperatures.
template <typename REAL>
inline REAL entropy_density_from_ne_mu(const REAL ne_mu)
{
    REAL s; // s = min(f, 1-f), always in [0, 0.5]
    if (ne_mu >= 0)
    {
        if (ne_mu > 700)
        {
            return 0;
        }
        s = 1 / (1 + std::exp(ne_mu));
    }
    else
    {
        if (ne_mu < -700)
        {
            return 0;
        }
        s = 1 / (1 + std::exp(-ne_mu));
    }
    if (s <= 0)
    {
        return 0;
    }
    return s * std::log(s) + (1 - s) * std::log1p(-s);
}
} // namespace

template <typename REAL>
Sto_Func<REAL>::Sto_Func()
{
    this->tem = Occupy::gaussian_parameter;
    this->mu  = static_cast<REAL>(0.0);
    this->t   = static_cast<REAL>(0.0);
    this->sigma = static_cast<REAL>(0.0);
    this->targ_e = static_cast<REAL>(0.0);
}

template <typename REAL>
void Sto_Func<REAL>::set_E_range(REAL* Emin_in, REAL* Emax_in)
{
    Emin = Emin_in;
    Emax = Emax_in;
}

template <typename REAL>
REAL Sto_Func<REAL>::root_fd(REAL e) const
{
    REAL e_mu = (e - mu) / this->tem;
    if (e_mu > 72) {
        return 0;
    } else {
        return 1 / sqrt(1 + exp(e_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::nroot_fd(REAL e) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL ne_mu = (e * DeltaE + Ebar - mu) / this->tem;
    if (ne_mu > 72) {
        return 0;
    } else {
        return 1 / sqrt(1 + exp(ne_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::fd(REAL e) const
{
    REAL e_mu = (e - mu) / this->tem;
    if (e_mu > 36) {
        return 0;
    } else {
        return 1 / (1 + exp(e_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::nfd(REAL e) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL ne_mu = (e * DeltaE + Ebar - mu) / this->tem;
    if (ne_mu > 36) {
        return 0;
    } else {
        return 1 / (1 + exp(ne_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::nxfd(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    REAL ne_mu = (e - mu) / this->tem;
    if (ne_mu > 36) {
        return 0;
    } else {
        return e / (1 + exp(ne_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::fdlnfd(REAL e) const
{
    const REAL e_mu = (e - mu) / this->tem;
    return entropy_density_from_ne_mu<REAL>(e_mu);
}

template <typename REAL>
REAL Sto_Func<REAL>::nfdlnfd(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    const REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem;
    return entropy_density_from_ne_mu<REAL>(ne_mu);
}

template <typename REAL>
REAL Sto_Func<REAL>::n_root_fdlnfd(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    const REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem;
    const REAL g = -entropy_density_from_ne_mu<REAL>(ne_mu);
    return (g > 0) ? std::sqrt(g) : 0;
}

template <typename REAL>
REAL Sto_Func<REAL>::nroot_mfd(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem;
    if (ne_mu < -72) {
        return 0;
    } else {
        return 1 / sqrt(1 + exp(-ne_mu));
}
}

template <typename REAL>
REAL Sto_Func<REAL>::ncos(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    return cos(e * t);
}

template <typename REAL>
REAL Sto_Func<REAL>::nsin(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    return sin(e * t);
}

template <typename REAL>
REAL Sto_Func<REAL>::n_sin(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    return -sin(e * t);
}

template <typename REAL>
REAL Sto_Func<REAL>::gauss(REAL e) const
{
    REAL a = pow((targ_e - e), 2) / 2.0 / pow(sigma, 2);
    if (a > 72) {
        return 0;
    } else {
        return exp(-a) / sqrt(TWOPI) / sigma;
}
}

template <typename REAL>
REAL Sto_Func<REAL>::ngauss(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    REAL a = pow((targ_e - e), 2) / 2.0 / pow(sigma, 2);
    if (a > 72) {
        return 0;
    } else {
        return exp(-a) / sqrt(TWOPI) / sigma;
}
}

template <typename REAL>
REAL Sto_Func<REAL>::nroot_gauss(REAL rawe) const
{
    REAL Ebar = (*Emin + *Emax) / 2;
    REAL DeltaE = (*Emax - *Emin) / 2;
    REAL e = rawe * DeltaE + Ebar;
    REAL a = pow((targ_e - e), 2) / 4.0 / pow(sigma, 2);
    if (a > 72) {
        return 0;
    } else {
        return exp(-a) / sqrt(sqrt(TWOPI) * sigma);
}
}

// we only have two examples: double and float.
template class Sto_Func<double>;
#ifdef __FLOAT_FFTW
template class Sto_Func<float>;
#endif
