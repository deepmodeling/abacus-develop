#include "sto_func.h"
#include "module_elecstate/occupy.h"
#define TWOPI 6.283185307179586

template<typename REAL>
Sto_Func<REAL>::Sto_Func()
{
    this->tem = Occupy::gaussian_parameter;
}

template<typename REAL>
REAL Sto_Func<REAL>:: root_fd(REAL e)
{
    REAL e_mu = (e - mu) / this->tem ;
    if(e_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(e_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: nroot_fd(REAL e)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL ne_mu = (e * DeltaE + Ebar - mu) / this->tem ;
    if(ne_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(ne_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: fd(REAL e)
{
    REAL e_mu = (e - mu) / this->tem ;
    if(e_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(e_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: nfd(REAL e)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL ne_mu = (e * DeltaE + Ebar - mu) / this->tem ;
    if(ne_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(ne_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: nxfd(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    REAL ne_mu = (e - mu) / this->tem ;
    if(ne_mu > 36)
        return 0;
    else
        return e / (1 + exp(ne_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: fdlnfd(REAL e)
{
    REAL e_mu = (e - mu) / this->tem;
    if (e_mu > 36)
        return 0;
    else if (e_mu < -36)
        return 0;
    else
    {
        REAL f = 1 / (1 + exp(e_mu));
        if (f == 0 || f == 1)
            return 0;
        else
            return (f * log(f) + (1.0 - f) * log(1.0 - f));
    }
}

template<typename REAL>
REAL Sto_Func<REAL>:: nfdlnfd(REAL rawe)
{
    REAL Ebar = (Emin + Emax) / 2;
    REAL DeltaE = (Emax - Emin) / 2;
    REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem;
    if (ne_mu > 36)
        return 0;
    else if (ne_mu < -36)
        return 0;
    else
    {
        REAL f = 1 / (1 + exp(ne_mu));
        if (f == 0 || f == 1)
            return 0;
        else
            return f * log(f) + (1 - f) * log(1 - f);
    }
}

template<typename REAL>
REAL Sto_Func<REAL>:: n_root_fdlnfd(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem ;
    if(ne_mu > 36)
        return 0;
    else if(ne_mu < -36)
        return 0;
    else
    {
        REAL f = 1 / (1 + exp(ne_mu));
        if (f == 0 || f == 1)
            return 0;
        else
            return sqrt(-f * log(f) - (1-f) * log(1-f));
    }
}

template<typename REAL>
REAL Sto_Func<REAL>::n_fd(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL ne_mu = (rawe * DeltaE + Ebar - mu) / this->tem ;
    if(ne_mu > 36)
        return 1;
    else
        return 1 - 1 / (1 + exp(ne_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: ncos(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    return cos(e * t);
}

template<typename REAL>
REAL Sto_Func<REAL>:: nsin(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    return sin(e * t);
}

template<typename REAL>
REAL Sto_Func<REAL>:: n_sin(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    return -sin(e * t);
}

template<typename REAL>
REAL Sto_Func<REAL>::gauss(REAL e)
{
    REAL a = pow((targ_e-e),2)/2.0/pow(sigma,2);
    if(a > 72)
        return 0;
    else
        return  exp(-a) /sqrt(TWOPI) / sigma ;
}

template<typename REAL>
REAL Sto_Func<REAL>::ngauss(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    REAL a = pow((targ_e-e),2)/2.0/pow(sigma,2);
    if(a > 72)
        return 0;
    else
        return  exp(-a) /sqrt(TWOPI) / sigma ;
}

template<typename REAL>
REAL Sto_Func<REAL>::nroot_gauss(REAL rawe)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL e = rawe * DeltaE + Ebar;
    REAL a = pow((targ_e-e),2)/4.0/pow(sigma,2);
    if(a > 72)
        return 0;
    else
        return  exp(-a) /sqrt(sqrt(TWOPI) * sigma) ;
}

//we only have two examples: double and float.
template class Sto_Func<double>;
#ifdef __ENABLE_FLOAT_FFTW
template class Sto_Func<float>;
#endif
