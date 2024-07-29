#include "broyden_mixing.h"

#include "module_base/lapack_connector.h"
#include "module_base/memory.h"
#include "module_base/module_container/base/third_party/blas.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
namespace Base_Mixing
{
template void Broyden_Mixing::tem_push_data(Mixing_Data& mdata,
                                            const double* data_in,
                                            const double* data_out,
                                            std::function<void(double*)> screen,
                                            std::function<void(double*, const double*, const double*)> mix,
                                            const bool& need_calcoef);
template void Broyden_Mixing::tem_push_data(
    Mixing_Data& mdata,
    const std::complex<double>* data_in,
    const std::complex<double>* data_out,
    std::function<void(std::complex<double>*)> screen,
    std::function<void(std::complex<double>*, const std::complex<double>*, const std::complex<double>*)> mix,
    const bool& need_calcoef);

template <class FPTYPE>
void Broyden_Mixing::tem_push_data(Mixing_Data& mdata,
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
        this->ndim_cal_dF = std::min(this->ndim_cal_dF + 1, this->mixing_ndim);
        start_dF = (this->start_dF + 1) % this->mixing_ndim;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
        for (int i = 0; i < length; ++i)
        {
            FP_F[i] = F_tmp[i];
            // dF{n} = F{n-1} - F{n} = -(F{n} - F{n-1})
            FP_dF[start_dF * length + i] -= FP_F[i];
        }
    }
};

template void Broyden_Mixing::tem_cal_coef(const Mixing_Data& mdata,
                                           std::function<double(double*, double*)> inner_product);
template void Broyden_Mixing::tem_cal_coef(
    const Mixing_Data& mdata,
    std::function<double(std::complex<double>*, std::complex<double>*)> inner_product);

template <class FPTYPE>
void Broyden_Mixing::tem_cal_coef(const Mixing_Data& mdata, std::function<double(FPTYPE*, FPTYPE*)> inner_product)
{
    ModuleBase::TITLE("Charge_Mixing", "Simplified_Broyden_mixing");
    ModuleBase::timer::tick("Charge", "Broyden_mixing");
    if (address != &mdata && address != nullptr)
        ModuleBase::WARNING_QUIT(
            "Broyden_mixing",
            "One Broyden_Mixing object can only bind one Mixing_Data object to calculate coefficients");
    const int length = mdata.length;
    FPTYPE* FP_dF = static_cast<FPTYPE*>(dF);
    FPTYPE* FP_F = static_cast<FPTYPE*>(F);
    if (ndim_cal_dF > 0)
    {
        ModuleBase::matrix beta_tmp(ndim_cal_dF, ndim_cal_dF);
        ModuleBase::matrix beta_tmp2(ndim_cal_dF, ndim_cal_dF);
        ModuleBase::matrix beta_tmp_diag(ndim_cal_dF, ndim_cal_dF);
        // beta(i, j) = <dF_i, dF_j>
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            FPTYPE* dFi = FP_dF + i * length;
            for (int j = i; j < ndim_cal_dF; ++j)
            {
                if (i != start_dF && j != start_dF)
                {
                    beta_tmp(i, j) = beta(i, j);
                }
                else
                {
                    FPTYPE* dFj = FP_dF + j * length;
                    beta(i, j) = beta_tmp(i, j) = inner_product(dFi, dFj);
                }
                if (j != i)
                {
                    beta_tmp(j, i) = beta_tmp(i, j);
                }
            }
        }
        // value tmp2
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            for (int j = 0; j < ndim_cal_dF; ++j)
            {
                beta_tmp2(i, j) = beta_tmp(i, j);
                beta_tmp_diag(i, j) = beta_tmp(i, j);
            }
        }
        // diag
        double* val = new double[ndim_cal_dF];
        diag(ndim_cal_dF, beta_tmp_diag.c, val);
        // output beta matrix
        std::cout << "----------------------beta matrix Begin---------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            for (int j = 0; j < ndim_cal_dF; ++j)
            {
                std::cout << beta_tmp(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "-----------------------------------------------------------" << std::endl;
        // output eigenvalues
        std::cout << "----------------------Eigenvalues Begin---------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            std::cout << val[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "-------------------------------------------------------------" << std::endl;

        double* work = new double[ndim_cal_dF];
        int* iwork = new int[ndim_cal_dF];
        char uu = 'U';
        int info;
        dsytrf_(&uu, &ndim_cal_dF, beta_tmp.c, &ndim_cal_dF, iwork, work, &ndim_cal_dF, &info);
        if (info != 0)
            ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when factorizing beta.");
        dsytri_(&uu, &ndim_cal_dF, beta_tmp.c, &ndim_cal_dF, iwork, work, &info);
        if (info != 0)
            ModuleBase::WARNING_QUIT("Charge_Mixing", "Error when DSYTRI beta.");
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            for (int j = i + 1; j < ndim_cal_dF; ++j)
            {
                beta_tmp(i, j) = beta_tmp(j, i);
            }
        }
        //output inverse beta matrix
        std::cout << "----------------------Inverse beta matrix Begin---------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            for (int j = 0; j < ndim_cal_dF; ++j)
            {
                std::cout << beta_tmp(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "---------------------------------------------------------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            FPTYPE* dFi = FP_dF + i * length;
            work[i] = inner_product(dFi, FP_F);
        }
        // gamma[i] = \sum_j beta_tmp(i,j) * work[j]
        std::vector<double> gamma(ndim_cal_dF);
        container::BlasConnector::gemv('N',
                                       ndim_cal_dF,
                                       ndim_cal_dF,
                                       1.0,
                                       beta_tmp.c,
                                       ndim_cal_dF,
                                       work,
                                       1,
                                       0.0,
                                       gamma.data(),
                                       1);
        // output gamma
        std::cout << "----------------------x=A^{-1}b Begin---------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            std::cout << gamma[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        coef[mdata.start] = 1 + gamma[dFindex_move(0)];
        for (int i = 1; i < ndim_cal_dF; ++i)
        {
            coef[mdata.index_move(-i)] = gamma[dFindex_move(-i)] - gamma[dFindex_move(-i + 1)];
        }
        coef[mdata.index_move(-ndim_cal_dF)] = -gamma[dFindex_move(-ndim_cal_dF + 1)];
        
        // try another way to calculate coef
        double* work2 = new double[ndim_cal_dF];
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            FPTYPE* dFi = FP_dF + i * length;
            work2[i] = inner_product(dFi, FP_F);
        }
        // use dsysv
        int m = 1;
        dsysv_(&uu, &ndim_cal_dF, &m, beta_tmp2.c, &ndim_cal_dF, iwork, work2, &ndim_cal_dF, work, &ndim_cal_dF, &info);

        //ouput work2
        std::cout << "----------------------Ax=b Begin---------------------" << std::endl;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            std::cout << work2[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;

        // calculate the diff between gamma and work2
        double diff = 0;
        for (int i = 0; i < ndim_cal_dF; ++i)
        {
            diff += std::abs(gamma[i] - work2[i]);
        }
        std::cout << "----------------------diff between gamma and work2---------------------" << std::endl;
        std::cout << diff << std::endl;
        std::cout << "-------------------------------------------------------------------------" << std::endl;

        delete[] work;
        delete[] iwork;
    }
    else
    {
        coef[0] = 1.0;
    }

    FPTYPE* dFnext = FP_dF + dFindex_move(1) * length;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(FPTYPE))
#endif
    for (int i = 0; i < length; ++i)
    {
        dFnext[i] = FP_F[i];
    }
    ModuleBase::timer::tick("Charge", "Broyden_mixing");
};

void Broyden_Mixing::diag(int n, double* A, double* val) {
    // diagonalize a real symmetric matrix A
    // A: input matrix (n x n); overwritten with eigenvectors
    // val: eigenvalues
    // work space query
    int lwork = -1;
    int info = 0;
    std::vector<double> work(1);
    dsyev_("V", "U", &n, A, &n, val, work.data(), &lwork, &info);
    lwork = work[0];
    work.resize(lwork);
    // actual computation
    dsyev_("V", "U", &n, A, &n, val, work.data(), &lwork, &info);
}

} // namespace Base_Mixing
