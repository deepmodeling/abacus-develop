#include "plain_mixing.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
namespace Base_Mixing
{
template void Plain_Mixing::tem_push_data(Mixing_Data& mdata,
                                          const double* data_in,
                                          const double* data_out,
                                          std::function<void(double*)> screen,
                                          std::function<void(double*, const double*, const double*)> mix,
                                          const bool& need_calcoef);
template void Plain_Mixing::tem_push_data(
    Mixing_Data& mdata,
    const std::complex<double>* data_in,
    const std::complex<double>* data_out,
    std::function<void(std::complex<double>*)> screen,
    std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)> mix,
    const bool& need_calcoef);

template <class FPTYPE>
void Plain_Mixing::tem_push_data(Mixing_Data& mdata,
                                 const FPTYPE* data_in,
                                 const FPTYPE* data_out,
                                 std::function<void(FPTYPE*)> screen,
                                 std::function<void(FPTYPE*, const FPTYPE*, const FPTYPE*)> mix,
                                 const bool& need_calcoef)
{
    const size_t length = mdata.length;
    std::vector<FPTYPE> F_tmp(length);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int i = 0; i < length; ++i)
    {
        F_tmp[i] = data_out[i] - data_in[i];
    }

    // get screened F
    if (screen != nullptr)
        screen(F_tmp.data());

    // container::Tensor data = data_in + mixing_beta * F;
    std::vector<FPTYPE> data(length);
    mix(data.data(), data_in, F_tmp.data());

    mdata.push(data.data());
};

template void Plain_Mixing::simple_mix(double* data_new,
                                       const double* data_in,
                                       const double* data_out,
                                       const int& length,
                                       std::function<void(double*)> screen);
template void Plain_Mixing::simple_mix(std::complex<double>* data_new,
                                       const std::complex<double>* data_in,
                                       const std::complex<double>* data_out,
                                       const int& length,
                                       std::function<void(std::complex<double>*)> screen);
template <class FPTYPE>
void Plain_Mixing::simple_mix(FPTYPE* data_new,
                              const FPTYPE* data_in,
                              const FPTYPE* data_out,
                              const int& length,
                              std::function<void(FPTYPE*)> screen)
{
    if (screen == nullptr)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int ig = 0; ig < length; ig++)
        {
            data_new[ig] = data_in[ig] + this->mixing_beta * (data_out[ig] - data_in[ig]);
        }
    }
    else
    {
        std::vector<FPTYPE> F_tmp(length);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int i = 0; i < length; ++i)
        {
            F_tmp[i] = data_out[i] - data_in[i];
        }
        screen(F_tmp.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int i = 0; i < length; ++i)
        {
            data_new[i] = data_in[i] + this->mixing_beta * F_tmp[i];
        }
    }
}

} // namespace Base_Mixing