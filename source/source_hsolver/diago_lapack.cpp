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
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
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
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
}

#ifdef __MPI
 template<>
    void DiagoLapack<double>::diag_pool(hamilt::MatrixBlock<double>& h_mat,
    hamilt::MatrixBlock<double>& s_mat,
    psi::Psi<double>& psi,
    Real* eigenvalue_in,
    MPI_Comm& comm)
{
    ModuleBase::TITLE("DiagoLapack", "diag_pool");
    assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);
    this->dsygvx_diag(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi);
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
}
    template<>
    void DiagoLapack<std::complex<double>>::diag_pool(hamilt::MatrixBlock<std::complex<double>>& h_mat,
    hamilt::MatrixBlock<std::complex<double>>& s_mat,
    psi::Psi<std::complex<double>>& psi,
    Real* eigenvalue_in,
    MPI_Comm& comm)
{
    ModuleBase::TITLE("DiagoLapack", "diag_pool");
    assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);
    std::vector<double> eigen(PARAM.globalv.nlocal, 0.0);
    this->zhegvx_diag(h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigen.data(), psi);
    const int inc = 1;
    BlasConnector::copy(PARAM.inp.nbands, eigen.data(), inc, eigenvalue_in, inc);
}
#endif

template <typename T>
std::pair<int, std::vector<int>> DiagoLapack<T>::dsygvx_once(const int ncol,
                                const int nrow,
                                const double* const h_mat,
                                const double* const s_mat,
                                double* const ekb,
                                psi::Psi<double>& wfc_2d) const
{
    ModuleBase::matrix h_tmp(ncol, nrow, false);
    memcpy(h_tmp.c, h_mat, sizeof(double) * ncol * nrow);
    ModuleBase::matrix s_tmp(ncol, nrow, false);
    memcpy(s_tmp.c, s_mat, sizeof(double) * ncol * nrow);

    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = PARAM.inp.nbands, one = 1;
    int M = 0, NZ = 0, lwork = -1, liwork = -1, info = 0;
    double vl = 0, vu = 0;
    const double abstol = 0, orfac = -1;
    std::vector<double> work(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(PARAM.globalv.nlocal, 0);
    std::vector<int> iclustr(2 * GlobalV::DSIZE);
    std::vector<double> gap(GlobalV::DSIZE);

    // LAPACK dsygvx signature:
    // (ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, VL, VU, IL, IU,
    //  ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO)
    int n = PARAM.globalv.nlocal;
    int lda = n, ldb = n, ldz = n;
    dsygvx_(&itype,
        &jobz,
        &range,
        &uplo,
        &n,
        h_tmp.c,
        &lda,
        s_tmp.c,
        &ldb,
        &vl,
        &vu,
        &il,
        &iu,
        &abstol,
        &M,
        ekb,
        wfc_2d.get_pointer(),
        &ldz,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);
    if (info) {
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + std::string(__FILE__) + " line "
                                 + std::to_string(__LINE__));
    }

    // Query returned optimal lwork in work[0]
    lwork = static_cast<int>(work[0]);
    work.resize(std::max(lwork, 3), 0);
    // LAPACK integer workspace: use conservative size (5*N)
    liwork = std::max(1, 5 * n);
    iwork.resize(liwork, 0);

    dsygvx_(&itype,
        &jobz,
        &range,
        &uplo,
        &n,
        h_tmp.c,
        &lda,
        s_tmp.c,
        &ldb,
        &vl,
        &vu,
        &il,
        &iu,
        &abstol,
        &M,
        ekb,
        wfc_2d.get_pointer(),
        &ldz,
        work.data(),
        &lwork,
        iwork.data(),
        ifail.data(),
        &info);
    //	GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

    if (info == 0) {
        return std::make_pair(info, std::vector<int>{});
    } else if (info < 0) {
        return std::make_pair(info, std::vector<int>{});
    } else if (info % 2) {
        return std::make_pair(info, ifail);
    } else if (info / 2 % 2) {
        return std::make_pair(info, iclustr);
    } else if (info / 4 % 2) {
        return std::make_pair(info, std::vector<int>{M, NZ});
    } else if (info / 16 % 2) {
        return std::make_pair(info, ifail);
    } else {
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + std::string(__FILE__) + " line "
                                 + std::to_string(__LINE__));
    }
}

template <typename T>
std::pair<int, std::vector<int>> DiagoLapack<T>::zhegvx_once(const int ncol,
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
    int M = 0, NZ = 0, lwork = -1, lrwork = -1, liwork = -1, info = 0;
    const double abstol = 0, orfac = -1;
    
    const double vl = 0, vu = 0;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<double> rwork(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(PARAM.globalv.nlocal, 0);
    std::vector<int> iclustr(2 * GlobalV::DSIZE);
    std::vector<double> gap(GlobalV::DSIZE);

    // LAPACK zhegvx signature:
    // (ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, VL, VU, IL, IU,
    //  ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO)
    int n = PARAM.globalv.nlocal;
    int lda = n, ldb = n, ldz = n;
    zhegvx_(&itype,
        &jobz,
        &range,
        &uplo,
        &n,
        h_tmp.c,
        &lda,
        s_tmp.c,
        &ldb,
        &vl,
        &vu,
        &il,
        &iu,
        &abstol,
        &M,
        ekb,
        wfc_2d.get_pointer(),
        &ldz,
        work.data(),
        &lwork,
        rwork.data(),
        iwork.data(),
        ifail.data(),
        &info);
    if (info) {
        throw std::runtime_error("info=" + ModuleBase::GlobalFunc::TO_STRING(info) + ". "
                                 + std::string(__FILE__) + " line "
                                 + std::to_string(__LINE__));
    }

    // Query returned optimal lwork in work[0]
    lwork = static_cast<int>(work[0].real());
    work.resize(std::max(lwork, 1), 0);
    // rwork: choose conservative size. Use baseline 7*N plus degeneracy margin
    lrwork = std::max(3, 7 * n + this->degeneracy_max * n);
    rwork.resize(lrwork, 0);
    // LAPACK integer workspace: use conservative size (5*N)
    liwork = std::max(1, 5 * n);
    iwork.resize(liwork, 0);

    zhegvx_(&itype,
        &jobz,
        &range,
        &uplo,
        &n,
        h_tmp.c,
        &lda,
        s_tmp.c,
        &ldb,
        &vl,
        &vu,
        &il,
        &iu,
        &abstol,
        &M,
        ekb,
        wfc_2d.get_pointer(),
        &ldz,
        work.data(),
        &lwork,
        rwork.data(),
        iwork.data(),
        ifail.data(),
        &info);
    //	GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

    if (info == 0) {
        return std::make_pair(info, std::vector<int>{});
    } else if (info < 0) {
        return std::make_pair(info, std::vector<int>{});
    } else if (info % 2) {
        return std::make_pair(info, ifail);
    } else if (info / 2 % 2) {
        return std::make_pair(info, iclustr);
    } else if (info / 4 % 2) {
        return std::make_pair(info, std::vector<int>{M, NZ});
    } else if (info / 16 % 2) {
        return std::make_pair(info, ifail);
    } else {
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + std::string(__FILE__) + " line "
                                 + std::to_string(__LINE__));
    }
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
        const std::pair<int, std::vector<int>> info_vec = dsygvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        post_processing(info_vec.first, info_vec.second);
        if (info_vec.first == 0) {
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
        const std::pair<int, std::vector<int>> info_vec = zhegvx_once(ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        post_processing(info_vec.first, info_vec.second);
        if (info_vec.first == 0) {
            break;
        }
    }
}

template <typename T>
void DiagoLapack<T>::post_processing(const int info, const std::vector<int>& vec)
{
    const std::string str_info = "info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n";
    const std::string str_FILE
        = std::string(__FILE__) + " line " + std::to_string(__LINE__) + ".\n";
    const std::string str_info_FILE = str_info + str_FILE;

    if (info == 0)
    {
        return;
    }
    else if (info < 0)
    {
        const int info_negative = -info;
        const std::string str_index
            = (info_negative > 100)
                  ? ModuleBase::GlobalFunc::TO_STRING(info_negative / 100) + "-th argument "
                        + ModuleBase::GlobalFunc::TO_STRING(info_negative % 100) + "-entry is illegal.\n"
                  : ModuleBase::GlobalFunc::TO_STRING(info_negative) + "-th argument is illegal.\n";
        throw std::runtime_error(str_info_FILE + str_index);
    }
    else if (info % 2)
    {
        std::string str_ifail = "ifail = ";
        for (const int i: vec) {
            str_ifail += ModuleBase::GlobalFunc::TO_STRING(i) + " ";
        }
        throw std::runtime_error(str_info_FILE + str_ifail);
    }
    else if (info / 2 % 2)
    {
        int degeneracy_need = 0;
        for (int irank = 0; irank < GlobalV::DSIZE; ++irank) {
            degeneracy_need = std::max(degeneracy_need, vec[2 * irank + 1] - vec[2 * irank]);
        }
        const std::string str_need = "degeneracy_need = " + ModuleBase::GlobalFunc::TO_STRING(degeneracy_need) + ".\n";
        const std::string str_saved
            = "degeneracy_saved = " + ModuleBase::GlobalFunc::TO_STRING(this->degeneracy_max) + ".\n";
        if (degeneracy_need <= this->degeneracy_max)
        {
            throw std::runtime_error(str_info_FILE + str_need + str_saved);
        }
        else
        {
            GlobalV::ofs_running << str_need << str_saved;
            this->degeneracy_max = degeneracy_need;
            return;
        }
    }
    else if (info / 4 % 2)
    {
        const std::string str_M = "M = " + ModuleBase::GlobalFunc::TO_STRING(vec[0]) + ".\n";
        const std::string str_NZ = "NZ = " + ModuleBase::GlobalFunc::TO_STRING(vec[1]) + ".\n";
        const std::string str_NBANDS
            = "PARAM.inp.nbands = " + ModuleBase::GlobalFunc::TO_STRING(PARAM.inp.nbands) + ".\n";
        throw std::runtime_error(str_info_FILE + str_M + str_NZ + str_NBANDS);
    }
    else if (info / 16 % 2)
    {
        const std::string str_npos = "not positive definite = " + ModuleBase::GlobalFunc::TO_STRING(vec[0]) + ".\n";
        throw std::runtime_error(str_info_FILE + str_npos);
    }
    else
    {
        throw std::runtime_error(str_info_FILE);
    }
}
} // namespace hsolver