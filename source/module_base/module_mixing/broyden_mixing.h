#ifndef BROYDEN_MIXING_H_
#define BROYDEN_MIXING_H_
#include "mixing.h"
#include "module_base/matrix.h"

namespace Base_Mixing
{
/**
 * @brief Simplified modified broyden_mixing method.
 *        Ref: D.D. Johnson PRB 38, 12807 (1988)
 *        Here the weight w0 of the error of the inverse Jacobian is set to 0 and the weight wn of
 *        the error of each previous iteration is set to same.
 * @note  Formula:
 *        F = n_out - n_in
 *        dF{i} = F_{i-1} - F{i}  //different from Ref
 *        dn_in{i} = n_in_{i-1} - n_in{i}  //different from Ref
 *        alpha{ij} = <dF{i}, dF{j}>
 *        beta{ij} = inv(alpha){ij}
 *        c{mk} = <dF{k}, F{m}>
 *        gamma{mn} = \sum_k c{mk} * beta{kn}
 *        n{m+1} = n_in{m} + mixing_beta*F{m} - \sum_n gamma{mn} * (dn_in{n} + mixing_beta*dF{n})
 *        mixing_data{i} = n_in{i} + mixing_beta*F{i}
 *        n{m+1} = \sum_i coef{i} * mixing_data{i}
 */
class Broyden_Mixing : public Mixing
{
  public:
    Broyden_Mixing(const int& mixing_ndim)
    {
        this->mixing_ndim = mixing_ndim;
        this->data_ndim = mixing_ndim + 1;
        this->coef = std::vector<double>(mixing_ndim + 1);
        this->beta = ModuleBase::matrix(mixing_ndim, mixing_ndim, true);
    }
    Broyden_Mixing(const int& mixing_ndim, const double& mixing_beta) : Broyden_Mixing(mixing_ndim)
    {
        this->mixing_beta = mixing_beta;
    }
    virtual ~Broyden_Mixing() override
    {
        if (F != nullptr)
            free(F);
        if (dF != nullptr)
            free(dF);
    };
    /**
     * @brief reset mixing
     *
     */
    virtual void reset() override
    {
        this->ndim_cal_dF = 0;
        this->start_dF = -1;
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

  private:
    // F = data_out - data_in
    void* F = nullptr;
    // dF = F_{n+1} - F_n
    void* dF = nullptr;
    // binded mixing_data
    Mixing_Data* address = nullptr;
    // beta_ij = <dF_i, dF_j>
    ModuleBase::matrix beta;
    // mixing_ndim = data_ndim - 1
    int mixing_ndim = -1;

    // start index for dF
    int start_dF = -1;
    // get the index of i-th dF vector
    int dFindex_move(const int& index)
    {
        return (start_dF + index + mixing_ndim) % mixing_ndim;
    }
    // the number of calculated dF
    int ndim_cal_dF = 0;
};
} // namespace Base_Mixing
#endif