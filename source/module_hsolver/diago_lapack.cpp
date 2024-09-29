// Refactored according to diago_scalapack
// This code will be futher refactored to remove the dependency of psi and hamilt
#include "module_parameter/parameter.h"

#include "diago_lapack.h"

#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{
template <>
void DiagoLapack<double>::diag(hamilt::Hamilt<double>* phm_in, psi::Psi<double>& psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoLapack", "diag");
    // Prepare H and S matrix
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);

    assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);

    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);

    // Diag
    this->dsygvx_diag(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi);
    // Copy result
    int size = eigen.size();
    for (int i = 0; i < size; i++)
    {
        eigenvalue_in[i] = eigen[i];
    }
}
template <>
void DiagoLapack<std::complex<double>>::diag(hamilt::Hamilt<std::complex<double>>* phm_in,
                                             psi::Psi<std::complex<double>>& psi,
                                             Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoLapack", "diag");
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);

    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);
    this->zhegvx_diag(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi);
    int size = eigen.size();
    for (int i = 0; i < size; i++)
    {
        eigenvalue_in[i] = eigen[i];
    }
}

template <typename T>
int DiagoLapack<T>::dsygvx_once(const int ncol,
                                const int nrow,
                                const double* const h_mat,
                                const double* const s_mat,
                                double* const ekb,
                                psi::Psi<double>& wfc_2d) const
{
    // Copy matrix to temp variables
    ModuleBase::matrix h_tmp(ncol, nrow, false);
    memcpy(h_tmp.c, h_mat, sizeof(double) * ncol * nrow);


    ModuleBase::matrix s_tmp(ncol, nrow, false);
    memcpy(s_tmp.c, s_mat, sizeof(double) * ncol * nrow);

    // Prepare caculate parameters
    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = PARAM.inp.nbands, one = 1;
    int M = 0, info = 0;
    double vl = 0, vu = 0;
    const double abstol = 0;

    int lwork = (ncol + 2) * ncol;

    std::vector<double> work(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(PARAM.globalv.nlocal, 0);

    // Original Lapack caculate, obelsete
    /*dsygvx_(&itype,
            &jobz,
            &range,
            &uplo,
            &PARAM.globalv.nlocal,
            h_tmp.c,
            &ncol,
            s_tmp.c,
            &ncol,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &M,
            ekb,
            wfc_2d.get_pointer(),
            &ncol,
            work.data(),
            &lwork,
            iwork.data(),
            ifail.data(),
            &info);

    // Throw error if it returns info
    if (info)
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
    //lwork = work[0];
    //work.resize(std::max(lwork, 3), 0);
    //iwork.resize(iwork[0], 0);

    dsygvx_(&itype,
            &jobz,
            &range,
            &uplo,
            &PARAM.globalv.nlocal,
            h_tmp.c,
            &PARAM.globalv.nlocal,
            s_tmp.c,
            &PARAM.globalv.nlocal,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &M,
            ekb,
            wfc_2d.get_pointer(),
            &ncol,
            work.data(),
            &lwork,
            iwork.data(),
            ifail.data(),
            &info);*/

    double *ev = new double[ncol * ncol];

    dsygv_(&itype, &jobz, &uplo, &PARAM.globalv.nlocal, h_tmp.c, &ncol, s_tmp.c, &ncol, ekb, ev, &lwork, &info);

    return info;
}

template <typename T>
int DiagoLapack<T>::zhegvx_once(const int ncol,
                                const int nrow,
                                const std::complex<double>* const h_mat,
                                const std::complex<double>* const s_mat,
                                double* const ekb,
                                psi::Psi<std::complex<double>>& wfc_2d) const
{
    ModuleBase::ComplexMatrix h_tmp(ncol, nrow, false);
    memcpy(h_tmp.c, h_mat, sizeof(std::complex<double>) * ncol * nrow);

    ModuleBase::ComplexMatrix s_tmp(ncol, nrow, false);
    memcpy(s_tmp.c, s_mat, sizeof(std::complex<double>) * ncol * nrow);

    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = PARAM.inp.nbands, one = 1;
    int M = 0, lrwork = -1, info = 0;
    const double abstol = 0;

    int lwork = (ncol + 2) * ncol;

    const double vl = 0, vu = 0;
    std::vector<std::complex<double>> work(1, 0);
    double *rwork = new double[3 * ncol - 2];
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(PARAM.globalv.nlocal, 0);

    // Original Lapack caculate, obelsete
    /*
    zhegvx_(&itype,
            &jobz,
            &range,
            &uplo,
            &PARAM.globalv.nlocal,
            h_tmp.c,
            &PARAM.globalv.nlocal,
            s_tmp.c,
            &PARAM.globalv.nlocal,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &M,
            ekb,
            wfc_2d.get_pointer(),
            &ncol,
            work.data(),
            &lwork,
            rwork.data(),
            iwork.data(),
            ifail.data(),
            &info);

    if (info)
        throw std::runtime_error("info=" + ModuleBase::GlobalFunc::TO_STRING(info) + ". "
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));

    //	GlobalV::ofs_running<<"lwork="<<work[0]<<"\t"<<"lrwork="<<rwork[0]<<"\t"<<"liwork="<<iwork[0]<<std::endl;

    //lwork = work[0].real();
    //work.resize(lwork, 0);
    //int maxlrwork = std::max(lrwork, 3);
    //rwork.resize(maxlrwork, 0);
    //iwork.resize(iwork[0], 0);

    zhegvx_(&itype,
            &jobz,
            &range,
            &uplo,
            &PARAM.globalv.nlocal,
            h_tmp.c,
            &PARAM.globalv.nlocal,
            s_tmp.c,
            &PARAM.globalv.nlocal,
            &vl,
            &vu,
            &il,
            &iu,
            &abstol,
            &M,
            ekb,
            wfc_2d.get_pointer(),
            &ncol,
            work.data(),
            &lwork,
            rwork.data(),
            iwork.data(),
            ifail.data(),
            &info);

    */

    std::complex<double> *ev = new std::complex<double>[ncol * ncol];

    zhegv_(&itype, &jobz, &uplo, &PARAM.globalv.nlocal, h_tmp.c, &ncol, s_tmp.c, &ncol, ekb, ev, &lwork, rwork, &info);

    return info;
}

template <typename T>
void DiagoLapack<T>::dsygvx_diag(const int ncol,
                                 const int nrow,
                                 const double* const h_mat,
                                 const double* const s_mat,
                                 double* const ekb,
                                 psi::Psi<double>& wfc_2d)
{
    while (true)
    {

        int info_result = dsygvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        if (info_result == 0) {
            break;
        }
    }
}

template <typename T>
void DiagoLapack<T>::zhegvx_diag(const int ncol,
                                 const int nrow,
                                 const std::complex<double>* const h_mat,
                                 const std::complex<double>* const s_mat,
                                 double* const ekb,
                                 psi::Psi<std::complex<double>>& wfc_2d)
{
    while (true)
    {
        int info_result = zhegvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        if (info_result == 0) {
            break;
        }
    }
}

template <typename T>
void DiagoLapack<T>::post_processing(const int info, const std::vector<int>& vec)
{
    const std::string str_info = "info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n";
    const std::string str_FILE
        = ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__) + ".\n";
    const std::string str_info_FILE = str_info + str_FILE;

    if (info == 0)
    {
        return;
    }
}
} // namespace hsolver