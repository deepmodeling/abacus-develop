#include "mixing.h"

#include "module_base/module_container/base/third_party/blas.h"
namespace Base_Mixing
{

void Mixing::push_data(Mixing_Data& mdata,
                       const double* data_in,
                       const double* data_out,
                       std::function<void(double*)> screen,
                       const bool& need_calcoef)
{
    const size_t length = mdata.length;
    this->push_data(
        mdata,
        data_in,
        data_out,
        screen,
        [this, length](double* out, const double* in, const double* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
            for (int i = 0; i < length; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
        },
        need_calcoef);
    return;
}

void Mixing::push_data(Mixing_Data& mdata,
                       const std::complex<double>* data_in,
                       const std::complex<double>* data_out,
                       std::function<void(std::complex<double>*)> screen,
                       const bool& need_calcoef)
{
    const size_t length = mdata.length;
    this->push_data(
        mdata,
        data_in,
        data_out,
        screen,
        [this, length](std::complex<double>* out, const std::complex<double>* in, const std::complex<double>* sres) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
            for (int i = 0; i < length; ++i)
            {
                out[i] = in[i] + this->mixing_beta * sres[i];
            }
        },
        need_calcoef);
    return;
}

void Mixing::mix_data(const Mixing_Data& mdata, double* data_mix)
{
    if (mdata.length <= 0)
        return;
    double* FP_data = static_cast<double*>(mdata.data);
    if (mdata.ndim_use == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 512)
#endif
        for (int i = 0; i < mdata.length; ++i)
            data_mix[i] = FP_data[i];
        return;
    }
    container::BlasConnector::gemv('N',
                                   mdata.length,
                                   mdata.ndim_use,
                                   1.0,
                                   FP_data,
                                   mdata.length,
                                   coef.data(),
                                   1,
                                   0.0,
                                   data_mix,
                                   1);
}
void Mixing::mix_data(const Mixing_Data& mdata, std::complex<double>* data_mix)
{
    if (mdata.length <= 0)
        return;
    std::complex<double>* FP_data = static_cast<std::complex<double>*>(mdata.data);
    if (mdata.ndim_use == 1)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 256)
#endif
        for (int i = 0; i < mdata.length; ++i)
            data_mix[i] = FP_data[i];
        return;
    }
    // conver coef to complex
    std::vector<std::complex<double>> coef_complex(coef.size());
    for (int i = 0; i < coef.size(); ++i)
        coef_complex[i] = coef[i];
    container::BlasConnector::gemv('N',
                                   mdata.length,
                                   mdata.ndim_use,
                                   1.0,
                                   FP_data,
                                   mdata.length,
                                   coef_complex.data(),
                                   1,
                                   0.0,
                                   data_mix,
                                   1);
}
} // namespace Base_Mixing