// Refactored according to diago_scalapack
// This code will be futher refactored to remove the dependency of psi and hamilt
#include "source_io/module_parameter/parameter.h"

#include "diago_lapack.h"

#include "source_base/global_variable.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_pw/module_pwdft/global.h"

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
void DiagoLapack<T>::dsygvx_once(const int ncol,
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
    const char jobz = 'V', uplo = 'U';
    const int itype = 1;

    std::vector<double> ev(ncol * ncol);

    LapackConnector::sygv(itype, jobz, uplo, PARAM.globalv.nlocal, h_tmp.c, ncol, s_tmp.c, ncol, ekb, ev.data());
}

template <typename T>
void DiagoLapack<T>::zhegvx_once(const int ncol,
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

    const char jobz = 'V', uplo = 'U';
    const int itype = 1;

    std::vector<std::complex<double>> ev(ncol * ncol);

    LapackConnector::hegv(LapackConnector::ColMajor, itype, jobz, uplo, PARAM.globalv.nlocal, h_tmp.c, ncol, s_tmp.c, ncol, ekb, ev.data());
}

template <typename T>
void DiagoLapack<T>::dsygvx_diag(const int ncol,
                                 const int nrow,
                                 const double* const h_mat,
                                 const double* const s_mat,
                                 double* const ekb,
                                 psi::Psi<double>& wfc_2d)
{
    dsygvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
}

template <typename T>
void DiagoLapack<T>::zhegvx_diag(const int ncol,
                                 const int nrow,
                                 const std::complex<double>* const h_mat,
                                 const std::complex<double>* const s_mat,
                                 double* const ekb,
                                 psi::Psi<std::complex<double>>& wfc_2d)
{
    zhegvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
}
} // namespace hsolver