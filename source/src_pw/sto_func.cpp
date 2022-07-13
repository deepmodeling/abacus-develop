#include "sto_func.h"
#include "occupy.h"

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
    if(ne_mu > 40)
        return 0;
    else
        return e / (1 + exp(ne_mu));
}

template<typename REAL>
REAL Sto_Func<REAL>:: fdlnfd(REAL e)
{
    REAL e_mu = (e - mu) / this->tem ;
    if(e_mu > 36)
        return 0;
    else if(e_mu < -36)
        return 0;
    else
    {
        REAL f = 1 / (1 + exp(e_mu));
        return (f * log(f) + (1.0-f) * log(1.0-f));
    }
}

template<typename REAL>
REAL Sto_Func<REAL>:: nfdlnfd(REAL e)
{
    REAL Ebar = (Emin + Emax)/2;
	REAL DeltaE = (Emax - Emin)/2;
    REAL ne_mu = (e * DeltaE + Ebar - mu) / this->tem ;
    if(ne_mu > 36)
        return 0;
    else if(ne_mu < -36)
        return 0;
    else
    {
        REAL f = 1 / (1 + exp(ne_mu));
        return f * log(f) + (1-f) * log(1-f);
    }
}

//we only have two examples: REAL and float.
template class Sto_Func<double>;
#ifdef __MIX_PRECISION
template class Sto_Func<float>;
#endif
