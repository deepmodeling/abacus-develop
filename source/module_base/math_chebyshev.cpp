#include "math_chebyshev.h"

#include "blas_connector.h"
#include "constants.h"
#include "global_function.h"
#include "module_base/module_container/ATen/kernels/blas.h"
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

template <typename REAL, typename Device>
Chebyshev<REAL, Device>::Chebyshev(const int norder_in) : fftw(2 * EXTEND * norder_in)
{
    this->norder = norder_in;
    norder2 = 2 * norder * EXTEND;
    if (this->norder < 1)
    {
        ModuleBase::WARNING_QUIT("Chebyshev", "The Chebyshev expansion order should be at least 1!");
    }
    coefr_cpu = new REAL[norder];
    coefc_cpu = new std::complex<REAL>[norder];
    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        resmem_var_op()(this->coef_real, norder);
        resmem_complex_op()(this->coef_complex, norder);
    }
    else
    {
        coef_real = coefr_cpu;
        coef_complex = coefc_cpu;
    }
    polytrace = new REAL[norder];

    // ndmin = ndmax = ndmax_in;
    getcoef_complex = false;
    getcoef_real = false;
}

template <typename REAL, typename Device>
Chebyshev<REAL, Device>::~Chebyshev()
{
    delete[] polytrace;
    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        delmem_var_op()(this->coef_real);
        delmem_complex_op()(this->coef_complex);
    }
    else
    {
        coef_real = nullptr;
        coef_complex = nullptr;
    }

    delete[] coefr_cpu;
    delete[] coefc_cpu;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::getpolyval(const REAL x, REAL* polyval, const int N)
{
    polyval[0] = 1;
    polyval[1] = x;
    for (int i = 2; i < N; ++i)
    {
        polyval[i] = 2 * x * polyval[i - 1] - polyval[i - 2];
    }
}
template <typename REAL, typename Device>
inline REAL Chebyshev<REAL, Device>::recurs(const REAL x, const REAL Tn, REAL const Tn_1)
{
    return 2 * x * Tn - Tn_1;
}

template <typename REAL, typename Device>
REAL Chebyshev<REAL, Device>::ddot_real(const std::complex<REAL>* psi_L,
                                        const std::complex<REAL>* psi_R,
                                        const int N,
                                        const int LDA,
                                        const int m)
{
    REAL result = 0;
    const base_device::DEVICE_CPU* cpu_ctx = {};
    if (N == LDA || m == 1)
    {
        int dim2 = 2 * N * m;
        REAL *pL, *pR;
        pL = (REAL*)psi_L;
        pR = (REAL*)psi_R;
        REAL* dot_device = nullptr;
        resmem_var_op()(dot_device, 1);
        container::kernels::blas_dot<REAL, ct_Device>()(dim2, pL, 1, pR, 1, dot_device);
        syncmem_var_d2h_op()(&result, dot_device, 1);
        delmem_var_op()(dot_device);
    }
    else
    {
        REAL *pL, *pR;
        pL = (REAL*)psi_L;
        pR = (REAL*)psi_R;
        REAL* dot_device = nullptr;
        resmem_var_op()(dot_device, 1);
        for (int i = 0; i < m; ++i)
        {
            int dim2 = 2 * N;
            container::kernels::blas_dot<REAL, ct_Device>()(dim2, pL, 1, pR, 1, dot_device);
            REAL result_temp = 0;
            syncmem_var_d2h_op()(&result_temp, dot_device, 1);
            result += result_temp;
            pL += 2 * LDA;
            pR += 2 * LDA;
        }
        delmem_var_op()(dot_device);
    }
    return result;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calcoef_real(std::function<REAL(REAL)> fun)
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
            coefr_cpu[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3;
        }
        else
        {
            coefr_cpu[i] = (cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3;
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
            coefr_cpu[i] += real(pcoef[i]) / norder2 * 1 / 3;
        }
        else
        {
            coefr_cpu[i] += real(pcoef[i]) / norder2 * 2 / 3;
        }
    }

    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        syncmem_var_h2d_op()(coef_real, coefr_cpu, norder);
    }

    getcoef_real = true;
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calcoef_complex(std::function<std::complex<REAL>(std::complex<REAL>)> fun)
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
            coefc_cpu[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coefc_cpu[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
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
            coefc_cpu[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coefc_cpu[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
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
            coefc_cpu[i].real(real(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coefc_cpu[i].real(real(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 2 / 3);
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
            coefc_cpu[i].imag(imag(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coefc_cpu[i].imag(imag(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }
    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        syncmem_complex_h2d_op()(coef_complex, coefc_cpu, norder);
    }

    getcoef_complex = true;
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calcoef_pair(std::function<REAL(REAL)> fun1, std::function<REAL(REAL)> fun2)
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
            coefc_cpu[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coefc_cpu[i].real((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
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
            coefc_cpu[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 2 / 3);
        }
        else
        {
            coefc_cpu[i].imag((cos(phi) * pcoef[i].real() + sin(phi) * pcoef[i].imag()) / norder2 * 4 / 3);
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
            coefc_cpu[i].real(real(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coefc_cpu[i].real(real(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 2 / 3);
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
            coefc_cpu[i].imag(imag(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 1 / 3);
        }
        else
        {
            coefc_cpu[i].imag(imag(coefc_cpu[i]) + real(pcoef[i]) / norder2 * 2 / 3);
        }
    }

    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        syncmem_complex_h2d_op()(coef_complex, coefc_cpu, norder);
    }

    getcoef_complex = true;
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calfinalvec_real(
    std::function<void(std::complex<REAL>*, std::complex<REAL>*, const int)> funA,
    std::complex<REAL>* wavein,
    std::complex<REAL>* waveout,
    const int N,
    const int LDA,
    const int m)
{
    if (!getcoef_real)
    {
        ModuleBase::WARNING_QUIT("Chebyshev<REAL>", "Please calculate coef_real first!");
    }

    std::complex<REAL>* arraynp1 = nullptr;
    std::complex<REAL>* arrayn = nullptr;
    std::complex<REAL>* arrayn_1 = nullptr;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1)
    {
        ndmxt = N * m;
    }
    else
    {
        ndmxt = LDA * m;
    }

    resmem_complex_op()(arraynp1, ndmxt);
    resmem_complex_op()(arrayn, ndmxt);
    resmem_complex_op()(arrayn_1, ndmxt);

    memcpy_complex_op()(arrayn_1, wavein, ndmxt);
    // ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

    funA(arrayn_1, arrayn, m);

    // 0- & 1-st order
    setmem_complex_op()(waveout, 0, ndmxt);
    std::complex<REAL> coef0 = std::complex<REAL>(coefr_cpu[0], 0);
    container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coef0, arrayn_1, 1, waveout, 1);
    std::complex<REAL> coef1 = std::complex<REAL>(coefr_cpu[1], 0);
    container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coef1, arrayn, 1, waveout, 1);
    // for (int i = 0; i < ndmxt; ++i)
    // {
    //     waveout[i] = coef_real[0] * arrayn_1[i] + coef_real[1] * arrayn[i];
    // }

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        std::complex<REAL> coefior = std::complex<REAL>(coefr_cpu[ior], 0);
        container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coefior, arraynp1, 1, waveout, 1);
        // for (int i = 0; i < ndmxt; ++i)
        // {
        //     waveout[i] += coef_real[ior] * arraynp1[i];
        // }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }
    delmem_complex_op()(arraynp1);
    delmem_complex_op()(arrayn);
    delmem_complex_op()(arrayn_1);
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calfinalvec_complex(
    std::function<void(std::complex<REAL>*, std::complex<REAL>*, const int)> funA,
    std::complex<REAL>* wavein,
    std::complex<REAL>* waveout,
    const int N,
    const int LDA,
    const int m)
{
    if (!getcoef_complex)
    {
        ModuleBase::WARNING_QUIT("Chebyshev", "Please calculate coef_complex first!");
    }

    std::complex<REAL>* arraynp1 = nullptr;
    std::complex<REAL>* arrayn = nullptr;
    std::complex<REAL>* arrayn_1 = nullptr;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1)
    {
        ndmxt = N * m;
    }
    else
    {
        ndmxt = LDA * m;
    }

    resmem_complex_op()(arraynp1, ndmxt);
    resmem_complex_op()(arrayn, ndmxt);
    resmem_complex_op()(arrayn_1, ndmxt);

    memcpy_complex_op()(arrayn_1, wavein, ndmxt);

    funA(arrayn_1, arrayn, m);

    // 0- & 1-st order
    setmem_complex_op()(waveout, 0, ndmxt);
    container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coefc_cpu[0], arrayn_1, 1, waveout, 1);
    container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coefc_cpu[1], arrayn, 1, waveout, 1);
    // for (int i = 0; i < ndmxt; ++i)
    // {
    //     waveout[i] = coef_complex[0] * arrayn_1[i] + coef_complex[1] * arrayn[i];
    // }

    // more than 1-st orders
    for (int ior = 2; ior < norder; ++ior)
    {
        recurs_complex(funA, arraynp1, arrayn, arrayn_1, N, LDA, m);
        container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(ndmxt, &coefc_cpu[ior], arraynp1, 1, waveout, 1);
        // for (int i = 0; i < ndmxt; ++i)
        // {
        //     waveout[i] += coef_complex[ior] * arraynp1[i];
        // }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }
    delmem_complex_op()(arraynp1);
    delmem_complex_op()(arrayn);
    delmem_complex_op()(arrayn_1);
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::calpolyvec_complex(
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
        memcpy_complex_op()(tmpout, tmpin, N);
        // ModuleBase::GlobalFunc::DCOPY(tmpin, tmpout, N);
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

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::tracepolyA(
    std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
    std::complex<REAL>* wavein,
    const int N,
    const int LDA,
    const int m)
{
    std::complex<REAL>* arraynp1 = nullptr;
    std::complex<REAL>* arrayn = nullptr;
    std::complex<REAL>* arrayn_1 = nullptr;
    assert(N >= 0 && LDA >= N);
    int ndmxt;
    if (m == 1)
    {
        ndmxt = N * m;
    }
    else
    {
        ndmxt = LDA * m;
    }

    resmem_complex_op()(arraynp1, ndmxt);
    resmem_complex_op()(arrayn, ndmxt);
    resmem_complex_op()(arrayn_1, ndmxt);

    memcpy_complex_op()(arrayn_1, wavein, ndmxt);
    // ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, ndmxt);

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

    delmem_complex_op()(arraynp1);
    delmem_complex_op()(arrayn);
    delmem_complex_op()(arrayn_1);
    return;
}

template <typename REAL, typename Device>
void Chebyshev<REAL, Device>::recurs_complex(
    std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
    std::complex<REAL>* arraynp1,
    std::complex<REAL>* arrayn,
    std::complex<REAL>* arrayn_1,
    const int N,
    const int LDA,
    const int m)
{
    funA(arrayn, arraynp1, m);
    const std::complex<REAL> two = 2.0;
    const std::complex<REAL> invone = -1.0;
    for (int ib = 0; ib < m; ++ib)
    {
        container::kernels::blas_scal<std::complex<REAL>, ct_Device>()(N, &two, arraynp1 + ib * LDA, 1);
        container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(N,
                                                                    &invone,
                                                                    arrayn_1 + ib * LDA,
                                                                    1,
                                                                    arraynp1 + ib * LDA,
                                                                    1);

        // for (int i = 0; i < N; ++i)
        // {
        //     arraynp1[i + ib * LDA] = REAL(2.0) * arraynp1[i + ib * LDA] - arrayn_1[i + ib * LDA];
        // }
    }
}

template <typename REAL, typename Device>
bool Chebyshev<REAL, Device>::checkconverge(
    std::function<void(std::complex<REAL>* in, std::complex<REAL>* out, const int)> funA,
    std::complex<REAL>* wavein,
    const int N,
    const int LDA,
    REAL& tmax,
    REAL& tmin,
    REAL stept)
{
    bool converge = true;
    std::complex<REAL>* arraynp1 = nullptr;
    std::complex<REAL>* arrayn = nullptr;
    std::complex<REAL>* arrayn_1 = nullptr;

    resmem_complex_op()(arraynp1, LDA);
    resmem_complex_op()(arrayn, LDA);
    resmem_complex_op()(arrayn_1, LDA);

    memcpy_complex_op()(arrayn_1, wavein, N);
    // ModuleBase::GlobalFunc::DCOPY(wavein, arrayn_1, N);

    if (tmin == tmax)
    {
        tmax += stept;
    }

    funA(arrayn_1, arrayn, 1);
    REAL sum1, sum2;
    REAL t;
    if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
    {
        sum1 = this->ddot_real(arrayn_1, arrayn_1, N);
        sum2 = this->ddot_real(arrayn_1, arrayn, N);
    }
    else
    {
#ifdef __MPI
        sum1 = ModuleBase::GlobalFunc::ddot_real(N, arrayn_1, arrayn_1);
        sum2 = ModuleBase::GlobalFunc::ddot_real(N, arrayn_1, arrayn);
#else
        sum1 = this->ddot_real(arrayn_1, arrayn_1, N);
        sum2 = this->ddot_real(arrayn_1, arrayn, N);
#endif
    }
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
        if (base_device::get_device_type<Device>(this->ctx) == base_device::GpuDevice)
        {
            sum1 = this->ddot_real(arrayn, arrayn, N);
            sum2 = this->ddot_real(arrayn, arraynp1, N);
        }
        else
        {
#ifdef __MPI
            sum1 = ModuleBase::GlobalFunc::ddot_real(N, arrayn, arrayn);
            sum2 = ModuleBase::GlobalFunc::ddot_real(N, arrayn, arraynp1);
#else
            sum1 = this->ddot_real(arrayn, arrayn, N);
            sum2 = this->ddot_real(arrayn, arraynp1, N);
#endif
        }
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
        std::complex<REAL> two = 2.0;
        std::complex<REAL> invone = -1.0;
        container::kernels::blas_scal<std::complex<REAL>, ct_Device>()(N, &two, arraynp1, 1);
        container::kernels::blas_axpy<std::complex<REAL>, ct_Device>()(N, &invone, arrayn_1, 1, arraynp1, 1);
        // for (int i = 0; i < N; ++i)
        // {
        //     arraynp1[i] = REAL(2.0) * arraynp1[i] - arrayn_1[i];
        // }
        std::complex<REAL>* tem = arrayn_1;
        arrayn_1 = arrayn;
        arrayn = arraynp1;
        arraynp1 = tem;
    }

    delmem_complex_op()(arraynp1);
    delmem_complex_op()(arrayn);
    delmem_complex_op()(arrayn_1);
    return converge;
}

// we only have two examples: double and float.
template class Chebyshev<double>;
#ifdef __ENABLE_FLOAT_FFTW
template class Chebyshev<float>;
#endif
#if ((defined __CUDA) || (defined __ROCM))
template class Chebyshev<double, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleBase
