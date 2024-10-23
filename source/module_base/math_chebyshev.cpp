#include "math_chebyshev.h"

#include "blas_connector.h"
#include "constants.h"
#include "global_function.h"
#include "tool_quit.h"

#include <cassert>

namespace ModuleBase
{

FFTW<double>::FFTW(const int norder2_in)
{
    ccoef = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * norder2_in);
    dcoef = (double*)fftw_malloc(sizeof(double) * norder2_in);
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
    ccoef = (fftwf_complex*)fftw_malloc(sizeof(fftwf_complex) * norder2_in);
    dcoef = (float*)fftw_malloc(sizeof(float) * norder2_in);
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

// A number to control the number of grids in C_n integration
#define EXTEND 16

template <typename REAL>
Chebyshev<REAL>::Chebyshev(const int norder_in) : fftw(2 * EXTEND * norder_in)
{
    this->norder = norder_in;
    norder2 = 2 * norder * EXTEND;
    if (this->norder < 1)
    {
        ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "The Chebyshev expansion order should be at least 1!");
    }
    polytrace = new REAL[norder];
    coef_real = new REAL[norder];
    coef_complex = new std::complex<REAL>[norder];

    // ndmin = ndmax = ndmax_in;

    getcoef_complex = false;
    getcoef_real = false;
}

template <typename REAL>
Chebyshev<REAL>::~Chebyshev()
{
    delete[] polytrace;
    delete[] coef_real;
    delete[] coef_complex;
}

template <typename REAL>
void Chebyshev<REAL>::getpolyval(const REAL x, REAL* polyval, const int N)
{
    polyval[0] = 1;
    polyval[1] = x;
    for (int i = 2; i < N; ++i)
    {
        polyval[i] = 2 * x * polyval[i - 1] - polyval[i - 2];
    }
}
template <typename REAL>
inline REAL Chebyshev<REAL>::recurs(const REAL x, const REAL Tn, REAL const Tn_1)
{
    return 2 * x * Tn - Tn_1;
}

template <typename REAL>
REAL Chebyshev<REAL>::ddot_real(const std::complex<REAL>* psi_L,
                                const std::complex<REAL>* psi_R,
                                const int N,
                                const int LDA,
                                const int m)
{
    REAL result = 0;
    if (N == LDA || m == 1)
    {
        int dim2 = 2 * N * m;
        REAL *pL, *pR;
        pL = (REAL*)psi_L;
        pR = (REAL*)psi_R;
        result = BlasConnector::dot(dim2, pL, 1, pR, 1);
    }
    else
    {
        REAL *pL, *pR;
        pL = (REAL*)psi_L;
        pR = (REAL*)psi_R;
        for (int i = 0; i < m; ++i)
        {
            int dim2 = 2 * N;
            result += BlasConnector::dot(dim2, pL, 1, pR, 1);
            pL += 2 * LDA;
            pR += 2 * LDA;
        }
    }
    return result;
}

template <typename REAL>
void Chebyshev<REAL>::calcoef_real(std::function<REAL(REAL)> fun)
{
    std::complex<REAL>* pcoef = (std::complex<REAL>*)this->fftw.ccoef;

    // three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun((REAL)cos((i + 0.5) * ModuleBase::TWO_PI / norder2));
    }

    // this->fftw.dcoef --FFT--> fftw.pcoef
    this->fftw.execute_fftw();

    for (int i = 0; i < norder; ++i)
    {
        REAL phi = i * ModuleBase::PI / norder2;
        if (i == 0)
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3;
        }
        else
        {
            coef_real[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3;
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun(cos(i * ModuleBase::TWO_PI / norder2));
    }

    // this->fftw.dcoef --FFT--> fftw.pcoef
    this->fftw.execute_fftw();

    for (int i = 0; i < norder; ++i)
    {
        if (i == 0)
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 1 / 3;
        }
        else
        {
            coef_real[i] += real(pcoef[i]) / norder2 * 2 / 3;
        }
    }

    getcoef_real = true;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::calcoef_complex(std::function<std::complex<REAL>(std::complex<REAL>)> fun)
{
    std::complex<REAL>* pcoef = (std::complex<REAL>*)this->fftw.ccoef;

    // three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun(cos((i + 0.5) * ModuleBase::TWO_PI / norder2)).real();
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        REAL phi = i * ModuleBase::PI / norder2;
        if (i == 0)
        {
            coef_complex[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coef_complex[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
        }
    }

    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun(cos((i + 0.5) * ModuleBase::TWO_PI / norder2)).imag();
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        REAL phi = i * ModuleBase::PI / norder2;
        if (i == 0)
        {
            coef_complex[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coef_complex[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun(cos(i * ModuleBase::TWO_PI / norder2)).real();
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        if (i == 0)
        {
            coef_complex[i].real(real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coef_complex[i].real(real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }

    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun(cos(i * ModuleBase::TWO_PI / norder2)).imag();
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        if (i == 0)
        {
            coef_complex[i].imag(imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coef_complex[i].imag(imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }

    getcoef_complex = true;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::calcoef_pair(std::function<REAL(REAL)> fun1, std::function<REAL(REAL)> fun2)
{
    std::complex<REAL>* pcoef = (std::complex<REAL>*)this->fftw.ccoef;

    // three point = 2/3 M + 1/3 T;
    //-----------------------------------------------
    //(M)iddle point integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun1(cos((i + 0.5) * ModuleBase::TWO_PI / norder2));
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        REAL phi = i * ModuleBase::PI / norder2;
        if (i == 0)
        {
            coef_complex[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coef_complex[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
        }
    }

    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun2(cos((i + 0.5) * ModuleBase::TWO_PI / norder2));
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        REAL phi = i * ModuleBase::PI / norder2;
        if (i == 0)
        {
            coef_complex[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coef_complex[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
        }
    }

    //-----------------------------------------------
    //(T)rapezoid integral method part
    //-----------------------------------------------
    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun1(cos(i * ModuleBase::TWO_PI / norder2));
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        if (i == 0)
        {
            coef_complex[i].real(real(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coef_complex[i].real(real(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }

    for (int i = 0; i < norder2; ++i)
    {
        this->fftw.dcoef[i] = fun2(cos(i * ModuleBase::TWO_PI / norder2));
    }
    this->fftw.execute_fftw();
    for (int i = 0; i < norder; ++i)
    {
        if (i == 0)
        {
            coef_complex[i].imag(imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coef_complex[i].imag(imag(coef_complex[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }

    getcoef_complex = true;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::calfinalvec_real(std::function<void(std::complex<REAL>*, std::complex<REAL>*, const int)> funA,
                                       std::complex<REAL>* wavein,
                                       std::complex<REAL>* waveout,
                                       const int N,
                                       const int LDA,
                                       const int m)
{
    if (!getcoef_real) {
        ModuleBase::WARNING_QUIT("Chebyshev<REAL>", "Please calculate coef_real first!");
}

    std::complex<REAL>* arraynp1;
    std::complex<REAL>* arrayn;
    std::complex<REAL>* arrayn_1;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1) {
        ndmxt = N * m;
    } else {
        ndmxt = LDA * m;
}

    arraynp1 = new std::complex<REAL>[ndmxt];
    arrayn = new std::complex<REAL>[ndmxt];
    arrayn_1 = new std::complex<REAL>[ndmxt];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

    funA(arrayn_1, arrayn, m);

    // 0- & 1-st order
    for (int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_real[0] * arrayn_1[i] + coef_real[1] * arrayn[i];
    }

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        for (int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_real[ior] * arraynp1[i];
        }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }
    delete[] arraynp1;
    delete[] arrayn;
    delete[] arrayn_1;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::calfinalvec_complex(std::function<void(std::complex<REAL>*, std::complex<REAL>*, const int)> funA,
                                          std::complex<REAL>* wavein,
                                          std::complex<REAL>* waveout,
                                          const int N,
                                          const int LDA,
                                          const int m)
{
    if (!getcoef_complex) {
        ModuleBase::WARNING_QUIT("Stochastic_Chebychev", "Please calculate coef_complex first!");
}

    std::complex<REAL>* arraynp1;
    std::complex<REAL>* arrayn;
    std::complex<REAL>* arrayn_1;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1) {
        ndmxt = N * m;
    } else {
        ndmxt = LDA * m;
}

    arraynp1 = new std::complex<REAL>[ndmxt];
    arrayn = new std::complex<REAL>[ndmxt];
    arrayn_1 = new std::complex<REAL>[ndmxt];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

    funA(arrayn_1, arrayn, m);

    // 0- & 1-st order
    for (int i = 0; i < ndmxt; ++i)
    {
        waveout[i] = coef_complex[0] * arrayn_1[i] + coef_complex[1] * arrayn[i];
    }

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        for (int i = 0; i < ndmxt; ++i)
        {
            waveout[i] += coef_complex[ior] * arraynp1[i];
        }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }
    delete[] arraynp1;
    delete[] arrayn;
    delete[] arrayn_1;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::calpolyvec_complex(
    std::function<void(std::complex<REAL>*, std::complex<REAL>*, const int)> funA,
    std::complex<REAL>* wavein,
    std::complex<REAL>* polywaveout,
    const int N,
    const int LDA,
    const int m)
{

    assert(N >= 0 && LDA >= N);
    const int ndmxt = LDA * m;

    std::complex<REAL>* arraynp1 = polywaveout + 2 * ndmxt;
    std::complex<REAL>* arrayn = polywaveout + ndmxt;
    std::complex<REAL>* arrayn_1 = polywaveout;

    std::complex<REAL>*tmpin = wavein, *tmpout = arrayn_1;
    for (int i = 0; i < m; ++i)
    {
        ModuleBase::GlobalFunc::DCOPY(tmpin, tmpout, N);
        tmpin += LDA;
        tmpout += LDA;
    }

    // 1-st order
    funA(arrayn_1, arrayn, m);

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        arrayn_1 += ndmxt;
        arrayn += ndmxt;
        arraynp1 += ndmxt;
    }
    return;
}

template <typename REAL>
void Chebyshev<REAL>::tracepolyA(std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
                                 std::complex<REAL>* wavein,
                                 const int N,
                                 const int LDA,
                                 const int m)
{
    std::complex<REAL>* arraynp1;
    std::complex<REAL>* arrayn;
    std::complex<REAL>* arrayn_1;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1) {
        ndmxt = N * m;
    } else {
        ndmxt = LDA * m;
}

    arraynp1 = new std::complex<REAL>[ndmxt];
    arrayn = new std::complex<REAL>[ndmxt];
    arrayn_1 = new std::complex<REAL>[ndmxt];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

    funA(arrayn_1, arrayn, m);

    polytrace[0] = this->ddot_real(wavein, wavein, N, LDA, m);
    polytrace[1] = this->ddot_real(wavein, arrayn, N, LDA, m);

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        polytrace[ior] = this->ddot_real(wavein, arraynp1, N, LDA, m);
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }

    delete[] arraynp1;
    delete[] arrayn;
    delete[] arrayn_1;
    return;
}

template <typename REAL>
void Chebyshev<REAL>::recurs_complex(
    std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
    std::complex<REAL>* arraynp1,
    std::complex<REAL>* arrayn,
    std::complex<REAL>* arrayn_1,
    const int N,
    const int LDA,
    const int m)
{
    funA(arrayn, arraynp1, m);
    for (int ib = 0; ib < m; ++ib)
    {
        for (int i = 0; i < N; ++i)
        {
            arraynp1[i + ib * LDA] = REAL(2.0) * arraynp1[i + ib * LDA] - arrayn_1[i + ib * LDA];
        }
    }
}

template <typename REAL>
bool Chebyshev<REAL>::checkconverge(
    std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
    std::complex<REAL>* wavein,
    const int N,
    const int LDA,
    REAL& tmax,
    REAL& tmin,
    REAL stept)
{
    bool converge = true;
    std::complex<REAL>* arraynp1;
    std::complex<REAL>* arrayn;
    std::complex<REAL>* arrayn_1;

    arraynp1 = new std::complex<REAL>[LDA];
    arrayn = new std::complex<REAL>[LDA];
    arrayn_1 = new std::complex<REAL>[LDA];

    ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, N);

    if (tmin == tmax)
    {
        tmax += stept;
    }

    funA(arrayn_1, arrayn, 1);
    REAL sum1, sum2;
    REAL t;
#ifdef __MPI
    sum1 = ModuleBase::GlobalFunc::ddot_real(N, arrayn_1, arrayn_1);
    sum2 = ModuleBase::GlobalFunc::ddot_real(N, arrayn_1, arrayn);
#else
    sum1 = this->ddot_real(arrayn_1, arrayn_1, N);
    sum2 = this->ddot_real(arrayn_1, arrayn, N);
#endif
    t = sum2 / sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
    if (t < tmin || tmin == 0)
    {
        converge = false;
        tmin = t - stept;
    }
    if (t > tmax)
    {
        converge = false;
        tmax = t + stept;
    }

    for (int ior = 2; ior < norder; ++ior)
    {
        funA(arrayn, arraynp1, 1);
#ifdef __MPI
        sum1 = ModuleBase::GlobalFunc::ddot_real(N, arrayn, arrayn);
        sum2 = ModuleBase::GlobalFunc::ddot_real(N, arrayn, arraynp1);
#else
        sum1 = this->ddot_real(arrayn, arrayn, N);
        sum2 = this->ddot_real(arrayn, arraynp1, N);
#endif
        t = sum2 / sum1 * (tmax - tmin) / 2 + (tmax + tmin) / 2;
        if (t < tmin)
        {
            converge = false;
            tmin = t - stept;
        }
        else if (t > tmax)
        {
            converge = false;
            tmax = t + stept;
        }
        for (int i = 0; i < N; ++i)
        {
            arraynp1[i] = REAL(2.0) * arraynp1[i] - arrayn_1[i];
        }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }

    delete[] arraynp1;
    delete[] arrayn;
    delete[] arrayn_1;
    return converge;
}

// we only have two examples: double and float.
template class Chebyshev<double>;
#ifdef __ENABLE_FLOAT_FFTW
template class Chebyshev<float>;
#endif

} // namespace ModuleBase
