#include "math_chebyshev.h"

#include "blas_connector.h"
#include "constants.h"
#include "global_function.h"
namespace ModuleBase
{

FFTW<double>::FFTW(const int norder2_in)
{
    ccoef = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * norder2_in);
    dcoef = (double *) fftw_malloc(sizeof(double) * norder2_in);
    coef_plan = fftw_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
}
FFTW<double>::~FFTW()
{
    fftw_destroy_plan(coef_plan);
	fftw_free(ccoef);
	fftw_free(dcoef);
}
void FFTW<double>::execute_fftw()
{
    fftw_execute(this->coef_plan);
}

#ifdef __ENABLE_FLOAT_FFTW
FFTW<float>::FFTW(const int norder2_in)
{
    ccoef = (fftwf_complex *) fftw_malloc(sizeof(fftwf_complex) * norder2_in);
    dcoef = (float *) fftw_malloc(sizeof(float) * norder2_in);
    coef_plan = fftwf_plan_dft_r2c_1d(norder2_in, dcoef, ccoef, FFTW_ESTIMATE);
}
FFTW<float>::~FFTW()
{
    fftwf_destroy_plan(coef_plan);
	fftw_free(ccoef);
	fftw_free(dcoef);
}
void FFTW<float>::execute_fftw()
{
    fftwf_execute(this->coef_plan);
}
#endif

//A number to control the number of grids in C_n integration
#define EXTEND 16

template<typename REAL>
Chebyshev<REAL>::Chebyshev(const int norder_in) : fftw(2 * EXTEND * norder_in)
{
    this->norder = norder_in;
    norder2 = 2 * norder * EXTEND;
    if(this->norder < 1)
    {
        ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "The Chebyshev expansion order should be at least 1!");
    }
    polytrace = new REAL [norder];
    coef_real = new REAL [norder];
    coef_complex = new std::complex<REAL> [norder];

    // ndmin = ndmax = ndmax_in;

    getcoef_complex = false;
    getcoef_real = false;
}

template<typename REAL>
Chebyshev<REAL>::~Chebyshev()
{
	delete [] polytrace;
	delete [] coef_real;
    delete [] coef_complex;
}

template<typename REAL>
void Chebyshev<REAL>::getpolyval(const REAL x, REAL* polyval, const int N)
{
    polyval[0] = 1;
    polyval[1] = x;
    for(int i = 2; i < N; ++i)
    {
        polyval[i] = 2 * x * polyval[i-1] - polyval[i-2];
    }
}
template<typename REAL>
inline REAL Chebyshev<REAL>::recurs(const REAL x, const REAL Tn, REAL const Tn_1)
{
    return 2*x*Tn-Tn_1;
}

template<typename REAL>
REAL Chebyshev<REAL>:: ddot_real(
    const std::complex<REAL>* psi_L,
    const std::complex<REAL>* psi_R,
    const int N, const int LDA, const int m)
{
    REAL result = 0;
    if(N == LDA || m==1)
    {
        int dim2=2 * N * m;
        REAL *pL,*pR;
        pL=(REAL *)psi_L;
        pR=(REAL *)psi_R;
        result=BlasConnector::dot(dim2,pL,1,pR,1);
    }
    else
    {
        REAL *pL,*pR;
        pL=(REAL *)psi_L;
        pR=(REAL *)psi_R;
        for(int i = 0 ; i < m ; ++i)
        {
            int dim2=2 * N;
            result +=  BlasConnector::dot(dim2,pL,1,pR,1);
            pL += 2 * LDA;
            pR += 2 * LDA;
        }
    }
    return result;
}

//we only have two examples: double and float.
template class Chebyshev<double>;
#ifdef __ENABLE_FLOAT_FFTW
template class Chebyshev<float>;
#endif

}  // namespace ModuleBase
