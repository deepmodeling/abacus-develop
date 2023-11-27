#ifndef PULAY_MIXING_H_
#define PULAY_MIXING_H_
#include "mixing.h"
#include "module_base/matrix.h"

namespace Base_Mixing
{
/**
 * @brief Pulay mixing method.
 *        Ref: Pulay P. Chemical Physics Letters, 1980, 73(2): 393-398.
 * @note  Formula:
 *        F = n_out - n_in
 *        alpha{ij} = <F{i}, F{j}>
 *        beta{ij} = inv(alpha){ij}
 *        coef{i} = \sum_j beta{ij} / \sum_k \sum_j beta{kj}
 *        mixing_data{i} = n_in{i} + mixing_beta*F{i}
 *        n{m+1} = \sum_n coef{i} * mixing_data{i}
 */
class Pulay_Mixing : public Mixing
{
  public:
    Pulay_Mixing(const int& mixing_ndim)
    {
        this->mixing_ndim = mixing_ndim;
        this->data_ndim = mixing_ndim;
        this->coef = std::vector<double>(mixing_ndim);
        this->beta = ModuleBase::matrix(mixing_ndim, mixing_ndim, true);
    }
    Pulay_Mixing(const int& mixing_ndim, const double& mixing_beta) : Pulay_Mixing(mixing_ndim)
    {
        this->mixing_beta = mixing_beta;
    }
    virtual ~Pulay_Mixing() override
    {
        if (F != nullptr)
            free(F);
    }
    /**
     * @brief reset mixing
     *
     */
    virtual void reset() override
    {
        this->start_F = 0;
        this->address = nullptr;
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
        tem_cal_coef(mdata, inner_product);
    }
    virtual void cal_coef(const Mixing_Data& mdata,
                          std::function<double(std::complex<double>*, std::complex<double>*)> inner_product) override
    {
        tem_cal_coef(mdata, inner_product);
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
     * @brief calculate coeficients for mixing
     *
     * @param mdata Mixing_Data
     * @param inner_product pointer to the inner dot function
     */
    template <class FPTYPE>
    void tem_cal_coef(const Mixing_Data& mdata, std::function<double(FPTYPE*, FPTYPE*)> inner_product);

    // F = data_out - data_in
    void* F = nullptr;
    // binded mixing_data
    Mixing_Data* address = nullptr;
    // beta_ij = <F_i, F_j>
    ModuleBase::matrix beta;
    // mixing_ndim = data_ndim - 1
    int mixing_ndim = -1;
    // start index for F
    int start_F = 0;
};
} // namespace Base_Mixing
#endif