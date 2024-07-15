#include "diago_scalapack_indep.h"

#include "module_base/global_variable.h"

#include <vector>

using namespace hsolver;

template <typename T>
DiagoScalapack<T>::DiagoScalapack(const std::vector<Real>& precondition_in,
                                                const int& nband_in,
                                                const int& nbasis_in,
                                                const double& diag_thr_in,
                                                const int& diag_nmax_in)
    : precondition(precondition_in), n_band(nband_in), dim(nbasis_in), nbase_x(nband_in * david_ndim_in),
      diag_thr(diag_thr_in), iter_nmax(diag_nmax_in), is_subspace(need_subspace_in), diag_comm(diag_comm_in)
{
    this->device = base_device::get_device_type<Device>(this->ctx);

    this->one = &this->cs.one;
    this->zero = &this->cs.zero;
    this->neg_one = &this->cs.neg_one;
}

template <typename T>
DiagoScalapack<T>::~DiagoScalapack()
{
}


template<>
void DiagoScalapack<double>::diag(double* h_mat, double* s_mat, const int ncol, const int nrow, const int* const desc, double* psi, Real* eigenvalue_in)
{
    // Desc is for h_mat
    ModuleBase::TITLE("DiagoScalapack", "diag");

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    this->pdsygvx_diag(desc, ncol, nrow, h_mat, s_mat, eigen.data(), psi);

    const int inc = 1;

    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
}


template<>
void DiagoScalapack<std::complex<double>>::diag(std::complex<double>* h_mat, std::complex<double>* s_mat, const int ncol, const int nrow, const int* const desc, std::complex<double>* psi, Real* eigenvalue_in)
{
    ModuleBase::TITLE("DiagoScalapack", "diag");

    std::vector<double> eigen(GlobalV::NLOCAL, 0.0);

    this->pzhegvx_diag(desc, ncol, nrow, h_mat, s_mat, eigen.data(), psi);

    const int inc = 1;

    BlasConnector::copy(GlobalV::NBANDS, eigen.data(), inc, eigenvalue_in, inc);
}


template <typename T>
std::pair<int, std::vector<int>> DiagoScalapack<T>::pdsygvx_once(const int* const desc,
                                                         const int ncol,
                                                         const int nrow,
                                                         const double *const h_mat,
                                                         const double *const s_mat,
                                                         double *const ekb,
                                                         double *wfc_2d) const
{

    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = GlobalV::NBANDS, one = 1;
    int M = 0, NZ = 0, lwork = -1, liwork = -1, info = 0;
    double vl = 0, vu = 0;
    const double abstol = 0, orfac = -1;
    std::vector<double> work(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(GlobalV::NLOCAL, 0);
    std::vector<int> iclustr(2 * GlobalV::DSIZE);
    std::vector<double> gap(GlobalV::DSIZE);

    pdsygvx_(&itype,
             &jobz,
             &range,
             &uplo,
             &GlobalV::NLOCAL,
             h_mat,
             &one,
             &one,
             desc,
             s_mat,
             &one,
             &one,
             desc,
             &vl,
             &vu,
             &il,
             &iu,
             &abstol,
             &M,
             &NZ,
             ekb,
             &orfac,
             wfc_2d,
             &one,
             &one,
             desc,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);
    if (info)
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));

    lwork = work[0];
    work.resize(std::max(lwork,3), 0);
    liwork = iwork[0];
    iwork.resize(liwork, 0);

    pdsygvx_(&itype,
             &jobz,
             &range,
             &uplo,
             &GlobalV::NLOCAL,
             h_mat,
             &one,
             &one,
             desc,
             s_mat,
             &one,
             &one,
             desc,
             &vl,
             &vu,
             &il,
             &iu,
             &abstol,
             &M,
             &NZ,
             ekb,
             &orfac,
             wfc_2d,
             &one,
             &one,
             desc,
             work.data(),
             &lwork,
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);

    if (info == 0)
        return std::make_pair(info, std::vector<int>{});
    else if (info < 0)
        return std::make_pair(info, std::vector<int>{});
    else if (info % 2)
        return std::make_pair(info, ifail);
    else if (info / 2 % 2)
        return std::make_pair(info, iclustr);
    else if (info / 4 % 2)
        return std::make_pair(info, std::vector<int>{M, NZ});
    else if (info / 16 % 2)
        return std::make_pair(info, ifail);
    else
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
}

template<typename T>
std::pair<int, std::vector<int>> DiagoScalapack<T>::pzhegvx_once(const int* const desc,
                                                         const int ncol,
                                                         const int nrow,
                                                         const std::complex<double> *const h_mat,
                                                         const std::complex<double> *const s_mat,
                                                         double *const ekb,
                                                         std::complex<double> *wfc_2d) const
{
    const char jobz = 'V', range = 'I', uplo = 'U';
    const int itype = 1, il = 1, iu = GlobalV::NBANDS, one = 1;
    int M = 0, NZ = 0, lwork = -1, lrwork = -1, liwork = -1, info = 0;
    const double abstol = 0, orfac = -1;
    //Note: pzhegvx_ has a bug
    //      We must give vl,vu a value, although we do not use range 'V'
    //      We must give rwork at least a memory of sizeof(double) * 3
    const double vl = 0, vu = 0;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<double> rwork(3, 0);
    std::vector<int> iwork(1, 0);
    std::vector<int> ifail(GlobalV::NLOCAL, 0);
    std::vector<int> iclustr(2 * GlobalV::DSIZE);
    std::vector<double> gap(GlobalV::DSIZE);

    pzhegvx_(&itype,
             &jobz,
             &range,
             &uplo,
             &GlobalV::NLOCAL,
             h_mat,
             &one,
             &one,
             desc,
             s_mat,
             &one,
             &one,
             desc,
             &vl,
             &vu,
             &il,
             &iu,
             &abstol,
             &M,
             &NZ,
             ekb,
             &orfac,
             wfc_2d,
             &one,
             &one,
             desc,
             work.data(),
             &lwork,
             rwork.data(),
             &lrwork,
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);
    if (info)
        throw std::runtime_error("info=" + ModuleBase::GlobalFunc::TO_STRING(info) + ". "
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));

    lwork = work[0].real();
    work.resize(lwork, 0);
    lrwork = rwork[0] + this->degeneracy_max * GlobalV::NLOCAL;
    int maxlrwork = std::max(lrwork,3);
    rwork.resize(maxlrwork, 0);
    liwork = iwork[0];
    iwork.resize(liwork, 0);

    pzhegvx_(&itype,
             &jobz,
             &range,
             &uplo,
             &GlobalV::NLOCAL,
             h_mat,
             &one,
             &one,
             desc,
             s_mat,
             &one,
             &one,
             desc,
             &vl,
             &vu,
             &il,
             &iu,
             &abstol,
             &M,
             &NZ,
             ekb,
             &orfac,
             wfc_2d,
             &one,
             &one,
             desc,
             work.data(),
             &lwork,
             rwork.data(),
             &lrwork,
             iwork.data(),
             &liwork,
             ifail.data(),
             iclustr.data(),
             gap.data(),
             &info);

    if (info == 0)
        return std::make_pair(info, std::vector<int>{});
    else if (info < 0)
        return std::make_pair(info, std::vector<int>{});
    else if (info % 2)
        return std::make_pair(info, ifail);
    else if (info / 2 % 2)
        return std::make_pair(info, iclustr);
    else if (info / 4 % 2)
        return std::make_pair(info, std::vector<int>{M, NZ});
    else if (info / 16 % 2)
        return std::make_pair(info, ifail);
    else
        throw std::runtime_error("info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n"
                                 + ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                 + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
}

template<typename T>
void DiagoScalapack<T>::pdsygvx_diag(const int* const desc,
                             const int ncol,
                             const int nrow,
                             const double *const h_mat,
                             const double *const s_mat,
                             double *const ekb,
                             double *wfc_2d)
{
    while (true)
    {
        const std::pair<int, std::vector<int>> info_vec = pdsygvx_once(desc, ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        post_processing(info_vec.first, info_vec.second);
        if (info_vec.first == 0)
            break;
    }
}

template<typename T>
void DiagoScalapack<T> ::pzhegvx_diag(const int* const desc,
                             const int ncol,
                             const int nrow,
                             const std::complex<double> *const h_mat,
                             const std::complex<double> *const s_mat,
                             double *const ekb,
                             std::complex<double> *wfc_2d)
{
    while (true)
    {
        const std::pair<int, std::vector<int>> info_vec = pzhegvx_once(desc, ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
        post_processing(info_vec.first, info_vec.second);
        if (info_vec.first == 0)
            break;
    }
}

    template<typename T>
    void DiagoScalapack<T>::post_processing(const int info, const std::vector<int>& vec)
{
    const std::string str_info = "info = " + ModuleBase::GlobalFunc::TO_STRING(info) + ".\n";
    const std::string str_FILE
        = ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line " + ModuleBase::GlobalFunc::TO_STRING(__LINE__) + ".\n";
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
        for (const int i: vec)
            str_ifail += ModuleBase::GlobalFunc::TO_STRING(i) + " ";
        throw std::runtime_error(str_info_FILE + str_ifail);
    }
    else if (info / 2 % 2)
    {
        int degeneracy_need = 0;
        for (int irank = 0; irank < GlobalV::DSIZE; ++irank)
            degeneracy_need = std::max(degeneracy_need, vec[2 * irank + 1] - vec[2 * irank]);
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
            = "GlobalV::NBANDS = " + ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS) + ".\n";
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