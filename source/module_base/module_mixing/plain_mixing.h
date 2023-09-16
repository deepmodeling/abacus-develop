#ifndef PLAIN_MIXING_H_
#define PLAIN_MIXING_H_
#include "mixing.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

namespace Base_Mixing
{
/**
 * @brief Plain mixing : rho_new = rho_in + mixing_beta * (rho_out - rho_in)
 *
 */
class Plain_Mixing : public Mixing
{
  public:
    Plain_Mixing()
    {
        this->mixing_ndim = 1;
        this->coef = std::vector<double>(1, 1.0);

    }
    virtual ~Plain_Mixing() override{};
    virtual void push_data(Mixing_Data& mdata,
                           const double* data_in,
                           const double* data_out,
                           std::function<void(double*)> screen,
                           const bool& need_calcoef) override
    {
        this->tem_push_data(mdata, data_in, data_out, screen, need_calcoef);
    };
    virtual void push_data(Mixing_Data& mdata,
                           const std::complex<double>* data_in,
                           const std::complex<double>* data_out,
                           std::function<void(std::complex<double>*)> screen,
                           const bool& need_calcoef) override
    {
        this->tem_push_data(mdata, data_in, data_out, screen, need_calcoef);
    };
    virtual void cal_coef(const Mixing_Data& mdata, std::function<double(double*, double*)> inner_dot) override
    {
        return;
    }
    virtual void cal_coef(const Mixing_Data& mdata,
                          std::function<double(std::complex<double>*, std::complex<double>*)> inner_dot) override
    {
        return;
    }

  private:
    template <class FPTYPE>
    void tem_push_data(Mixing_Data& mdata,
                       const FPTYPE* data_in,
                       const FPTYPE* data_out,
                       std::function<void(FPTYPE*)> screen,
                       const bool& need_calcoef)
    {
        const size_t length = mdata.length;
        std::vector<FPTYPE> F_tmp(length);

#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
#endif
        for (int i = 0; i < length; ++i)
        {
            data[i] = data_in[i] + this->mixing_beta * F_tmp[i];
        }

        mdata.push(data.data());
    };
};
} // namespace Base_Mixing
#endif