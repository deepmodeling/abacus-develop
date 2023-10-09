#ifndef PULAY_MIXING_H_
#define PULAY_MIXING_H_
#include "mixing.h"
#include "module_base/lapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

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
    Pulay_Mixing(const int& mixing_ndim, const double& mixing_beta)
    {
        this->mixing_ndim = mixing_ndim;
        this->data_ndim = mixing_ndim;
        this->mixing_beta = mixing_beta;
        this->coef = std::vector<double>(mixing_ndim);
        this->beta = ModuleBase::matrix(mixing_ndim, mixing_ndim, true);
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
     * @param need_calcoef whether need to calculate the coef
     *
     */
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
     * @param need_calcoef whether need to calculate the coef
     *
     */
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
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
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
                "Pulay_Mixing",
                "One Pulay_Mixing object can only bind one Mixing_Data object to calculate coefficients");

        FPTYPE* FP_F = static_cast<FPTYPE*>(F);
        if (mdata.ndim_use == 1)
        {
            address = &mdata;
            // allocate
            if (F != nullptr)
                free(F);
            F = malloc(sizeof(FPTYPE) * length * mixing_ndim);
            FP_F = static_cast<FPTYPE*>(F);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
            for (int i = 0; i < length; ++i)
            {
                FP_F[i] = F_tmp[i];
            }
        }
        else
        {
            start_F = (this->start_F + 1) % this->mixing_ndim;
            FPTYPE* FP_startF = FP_F + start_F * length;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
            for (int i = 0; i < length; ++i)
            {
                FP_startF[i] = F_tmp[i];
            }
        }
    };

    /**
     * @brief calculate coeficients for mixing
     *
     * @param mdata Mixing_Data
     * @param inner_product pointer to the inner dot function
     */
    template <class FPTYPE>
    void tem_cal_coef(const Mixing_Data& mdata, std::function<double(FPTYPE*, FPTYPE*)> inner_product)
    {
        ModuleBase::TITLE("Charge_Mixing", "Pulay_mixing");
        ModuleBase::timer::tick("Charge", "Pulay_mixing");
        if (address != &mdata && address != nullptr)
            ModuleBase::WARNING_QUIT(
                "Pulay_mixing",
                "One Pulay_Mixing object can only bind one Mixing_Data object to calculate coefficients");
        const int length = mdata.length;
        FPTYPE* FP_F = static_cast<FPTYPE*>(F);

        if (mdata.ndim_use > 1)
        {
            const int ndim_use = mdata.ndim_use;
            ModuleBase::matrix beta_tmp(ndim_use, ndim_use);
            // beta(i, j) = <F_i, F_j>
            for (int i = 0; i < ndim_use; ++i)
            {
                FPTYPE* Fi = FP_F + i * length;
                for (int j = i; j < ndim_use; ++j)
                {
                    if (i != start_F && j != start_F)
                    {
                        beta_tmp(i, j) = beta(i, j);
                    }
                    else
                    {
                        FPTYPE* Fj = FP_F + j * length;
                        beta(i, j) = beta_tmp(i, j) = inner_product(Fi, Fj);
                    }
                    if (j != i)
                    {
                        beta(j, i) = beta_tmp(j, i) = beta_tmp(i, j);
                    }
                }
            }

            double* work = new double[ndim_use];
            int* iwork = new int[ndim_use];
            char uu = 'U';
            int info;
            dsytrf_(&uu, &ndim_use, beta_tmp.c, &ndim_use, iwork, work, &ndim_use, &info);
            if (info != 0)
                ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when factorizing beta.");
            dsytri_(&uu, &ndim_use, beta_tmp.c, &ndim_use, iwork, work, &info);
            if (info != 0)
                ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when DSYTRI beta.");
            for (int i = 0; i < ndim_use; ++i)
            {
                for (int j = i + 1; j < ndim_use; ++j)
                {
                    beta_tmp(i, j) = beta_tmp(j, i);
                }
            }

            // coef{i} = \sum_j beta{ij} / \sum_k \sum_j beta{kj}
            double sum_beta = 0.;
            for (int i = 0; i < ndim_use; ++i)
            {
                for (int j = 0; j < ndim_use; ++j)
                {
                    sum_beta += beta_tmp(j, i);
                }
            }
            for (int i = 0; i < ndim_use; ++i)
            {
                coef[i] = 0.;
                for (int j = 0; j < ndim_use; ++j)
                {
                    coef[i] += beta_tmp(i, j);
                }
                coef[i] /= sum_beta;
            }
            delete[] work;
            delete[] iwork;
        }
        else
        {
            beta(0, 0) = inner_product(FP_F, FP_F);
            coef[0] = 1.0;
        }

        ModuleBase::timer::tick("Charge", "Pulay_mixing");
    };

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