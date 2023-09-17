#ifndef BROYDEN_MIXING_H_
#define BROYDEN_MIXING_H_
#include "mixing.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_base/global_function.h"

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
    Broyden_Mixing(const int& mixing_ndim, const double& mixing_beta)
    {
        this->mixing_ndim = mixing_ndim;
        this->mixing_beta = mixing_beta;
        this->coef = std::vector<double>(mixing_ndim);
        this->beta = ModuleBase::matrix(mixing_ndim, mixing_ndim, true);
    }
    virtual ~Broyden_Mixing() override
    {
        if (F != nullptr)
            free(F);
        if (dF != nullptr)
            free(dF);
    };
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
        tem_cal_coef(mdata, inner_dot);
    }
    virtual void cal_coef(const Mixing_Data& mdata,
                          std::function<double(std::complex<double>*, std::complex<double>*)> inner_dot) override
    {
        tem_cal_coef(mdata, inner_dot);
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

        if (!need_calcoef)
            return;

        if (address != &mdata && address != nullptr)
            ModuleBase::WARNING_QUIT(
                "Broyden_Mixing",
                "One Broyden_Mixing object can only bind one Mixing_Data object to calculate coefficients");

        FPTYPE* FP_dF = static_cast<FPTYPE*>(dF);
        FPTYPE* FP_F = static_cast<FPTYPE*>(F);
        if (mdata.ndim_use == 1)
        {
            address = &mdata;
            // allocate
            if (F != nullptr)
                free(F);
            F = malloc(sizeof(FPTYPE) * length);
            FP_F = static_cast<FPTYPE*>(F);
            if (dF != nullptr)
                free(dF);
            dF = malloc(sizeof(FPTYPE) * length * mixing_ndim);
            FP_dF = static_cast<FPTYPE*>(dF);

            ModuleBase::GlobalFunc::DCOPY(F_tmp.data(), FP_F, length);
        }
        else
        {
            const int previous = mdata.index_move(-1);
            ModuleBase::GlobalFunc::DCOPY(F_tmp.data(), FP_F, length);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
#endif
            for (int i = 0; i < length; ++i)
            {
                // dF{n} = F{n-1} - F{n} = -(F{n} - F{n-1})
                FP_dF[previous * length + i] -= FP_F[i];
            }
        }
    };

    template <class FPTYPE>
    void tem_cal_coef(const Mixing_Data& mdata, std::function<double(FPTYPE*, FPTYPE*)> inner_dot)
    {
        ModuleBase::TITLE("Charge_Mixing", "Simplified_Broyden_mixing");
        ModuleBase::timer::tick("Charge", "Broyden_mixing");
        if (address != &mdata && address != nullptr)
            ModuleBase::WARNING_QUIT(
                "Mixing",
                "One Mixing object can only bind one Mixing_Data object to calculate coefficients");
        const int length = mdata.length;
        const int start = mdata.start;
        FPTYPE* FP_dF = static_cast<FPTYPE*>(dF);
        FPTYPE* FP_F = static_cast<FPTYPE*>(F);
        if (mdata.ndim_previous > 0)
        {
            const int ndim_previous = mdata.ndim_previous;
            const int previous = mdata.index_move(-1);
            ModuleBase::matrix beta_tmp(ndim_previous, ndim_previous);
            // beta(i, j) = <dF_i, dF_j>
            for (int i = 0; i < ndim_previous; ++i)
            {
                FPTYPE* dFi = FP_dF + i * length;
                for (int j = i; j < ndim_previous; ++j)
                {
                    if (i != previous && j != previous)
                    {
                        beta_tmp(i, j) = beta(i, j);
                    }
                    FPTYPE* dFj = FP_dF + j * length;
                    beta(i, j) = beta_tmp(i, j) = inner_dot(dFi, dFj);
                    if (j != i)
                    {
                        beta(j, i) = beta_tmp(j, i) = beta_tmp(i, j);
                    }
                }
            }
            double* work = new double[ndim_previous];
            int* iwork = new int[ndim_previous];
            char uu = 'U';
            int info;
            dsytrf_(&uu, &ndim_previous, beta_tmp.c, &ndim_previous, iwork, work, &ndim_previous, &info);
            if (info != 0)
                ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when factorizing beta.");
            dsytri_(&uu, &ndim_previous, beta_tmp.c, &ndim_previous, iwork, work, &info);
            if (info != 0)
                ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when DSYTRI beta.");
            for (int i = 0; i < ndim_previous; ++i)
            {
                for (int j = i + 1; j < ndim_previous; ++j)
                {
                    beta_tmp(i, j) = beta_tmp(j, i);
                }
            }
            for (int i = 0; i < ndim_previous; ++i)
            {
                FPTYPE* dFi = FP_dF + i * length;
                work[i] = inner_dot(dFi, FP_F);
            }
            // gamma[i] = \sum_j beta_tmp(i,j) * work[j]
            std::vector<double> gamma(ndim_previous);
            container::BlasConnector::gemv('N',
                                           ndim_previous,
                                           ndim_previous,
                                           1.0,
                                           beta_tmp.c,
                                           ndim_previous,
                                           work,
                                           1,
                                           0.0,
                                           gamma.data(),
                                           1);

            coef[start] = 1 + gamma[mdata.index_move(-1)];
            for (int i = 1; i < ndim_previous - 1; ++i)
            {
                coef[mdata.index_move(-i)] = gamma[mdata.index_move(-i - 1)] - gamma[mdata.index_move(-i)];
            }
            coef[mdata.index_move(-ndim_previous)] = -gamma[mdata.index_move(-ndim_previous)];

            delete[] work;
            delete[] iwork;
        }
        else
        {
            coef[0] = 1.0;
        }

        FPTYPE* dFstart = FP_dF + start * length;
        ModuleBase::GlobalFunc::DCOPY(FP_F, dFstart, length);
        ModuleBase::timer::tick("Charge", "Broyden_mixing");
    };

  private:
    // F = data_out - data_in
    void* F = nullptr;
    // dF = F_{n+1} - F_n
    void* dF = nullptr;
    // binded mixing_data
    Mixing_Data* address = nullptr;
    // beta_ij = <dF_i, dF_j>
    ModuleBase::matrix beta;
};
} // namespace Base_Mixing
#endif