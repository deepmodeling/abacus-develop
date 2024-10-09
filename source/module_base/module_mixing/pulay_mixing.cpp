#include "pulay_mixing.h"

#include "module_base/lapack_connector.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
namespace Base_Mixing
{
template void Pulay_Mixing::tem_push_data(Mixing_Data& mdata,
                                          const double* data_in,
                                          const double* data_out,
                                          std::function<void(double*)> screen,
                                          std::function<void(double*, const double*, const double*)> mix,
                                          const bool& need_calcoef);
template void Pulay_Mixing::tem_push_data(
    Mixing_Data& mdata,
    const std::complex<double>* data_in,
    const std::complex<double>* data_out,
    std::function<void(std::complex<double>*)> screen,
    std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)> mix,
    const bool& need_calcoef);

template <class FPTYPE>
void Pulay_Mixing::tem_push_data(Mixing_Data& mdata,
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

template void Pulay_Mixing::tem_cal_coef(const Mixing_Data& mdata,
                                         std::function<double(double*, double*)> inner_product);
template void Pulay_Mixing::tem_cal_coef(
    const Mixing_Data& mdata,
    std::function<double(std::complex<double>*, std::complex<double>*)> inner_product);

template <class FPTYPE>
void Pulay_Mixing::tem_cal_coef(const Mixing_Data& mdata, std::function<double(FPTYPE*, FPTYPE*)> inner_product)
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
} // namespace Base_Mixing
