#ifndef PLAIN_MIXING_H_
#define PLAIN_MIXING_H_
#include "mixing.h"

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
        this->coef = std::vector<double>(1, 1.0);
    }
    Plain_Mixing(const double& mixing_beta) : Plain_Mixing()
    {
        this->mixing_beta = mixing_beta;
    }
    virtual ~Plain_Mixing() override{};

    /**
     * @brief reset mixing
     *
     */
    virtual void reset() override
    {
        return;
    }

    /**
     * @brief
     *
     * @param mdata store information of this iterative step
     * @param data_in x_in
     * @param data_out x_out = f(x_in)
     * @param screen pointer to the screen function for Ker-Ker mixing
     * @param mix (double* out, const double* in, const double* residual) calculate 'out' with 'in' and residual
     * @param need_calcoef whether need to calculate the coef
     *
     */
    virtual void push_data(Mixing_Data& mdata,
                           const double* data_in,
                           const double* data_out,
                           std::function<void(double*)> screen,
                           std::function<void(double*, const double*, const double*)> mix,
                           const bool& need_calcoef) override
    {
        this->tem_push_data(mdata, data_in, data_out, screen, mix, need_calcoef);
    };
    virtual void push_data(
        Mixing_Data& mdata,
        const std::complex<double>* data_in,
        const std::complex<double>* data_out,
        std::function<void(std::complex<double>*)> screen,
        std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)> mix,
        const bool& need_calcoef) override
    {
        this->tem_push_data(mdata, data_in, data_out, screen, mix, need_calcoef);
    };

    /**
     * @brief calculate coeficients for mixing
     *
     * @param mdata Mixing_Data
     * @param inner_product pointer to the inner dot function
     */
    virtual void cal_coef(const Mixing_Data& mdata, std::function<double(double*, double*)> inner_product) override
    {
        return;
    }
    virtual void cal_coef(const Mixing_Data& mdata,
                          std::function<double(std::complex<double>*, std::complex<double>*)> inner_product) override
    {
        return;
    }

    /**
     * @brief Simple plain mixing
     *        data_new = data_in + mixing_beta * (data_out - data_in)
     *
     * @param data_new can be the same as data_in or data_out
     */
    void plain_mix(double* data_new,
                  const double* data_in,
                  const double* data_out,
                  const int& length,
                  std::function<void(double*)> screen)
    {
        this->simple_mix(data_new, data_in, data_out, length, screen);
    }
    void plain_mix(std::complex<double>* data_new,
                  const std::complex<double>* data_in,
                  const std::complex<double>* data_out,
                  const int& length,
                  std::function<void(std::complex<double>*)> screen)
    {
        this->simple_mix(data_new, data_in, data_out, length, screen);
    }

  private:
    /**
     * @brief
     *
     * @param mdata store information of this iterative step
     * @param data_in x_in
     * @param data_out x_out = f(x_in)
     * @param screen pointer to the screen function for Ker-Ker mixing
     * @param mix (double* out, const double* in, const double* residual) calculate 'out' with 'in' and residual
     * @param need_calcoef whether need to calculate the coef
     *
     */
    template <class FPTYPE>
    void tem_push_data(Mixing_Data& mdata,
                       const FPTYPE* data_in,
                       const FPTYPE* data_out,
                       std::function<void(FPTYPE*)> screen,
                       std::function<void(FPTYPE*, const FPTYPE*, const FPTYPE*)> mix,
                       const bool& need_calcoef);

    /**
     * @brief Simple plain mixing
     *        data_new = data_in + mixing_beta * (data_out - data_in)
     *
     * @param data_new can be the same as data_in or data_out
     */
    template <class FPTYPE>
    void simple_mix(FPTYPE* data_new,
                    const FPTYPE* data_in,
                    const FPTYPE* data_out,
                    const int& length,
                    std::function<void(FPTYPE*)> screen);
};
} // namespace Base_Mixing
#endif