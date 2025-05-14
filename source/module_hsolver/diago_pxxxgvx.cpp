#include "diago_pxxxgvx.h"

#include "module_base/blacs_connector.h"
#include "module_base/scalapack_connector.h"

#include <complex>
#include <cstring>
#include <iostream>
#include <vector>

namespace hsolver
{

#ifdef __MPI

/**
 * @file diago_pxxxgvx.cpp
 * @brief This file contains functions for performing parallel diagonalization using Scalapack.
 *
 * The functions in this file are designed to handle different data types (double, float, std::complex<double>,
 * std::complex<float>) and perform diagonalization using Scalapack routines (pdsygvx_, pzhegvx_, pssygvx_, pchegvx_).
 *
 * The main functions are:
 * - pxxxgvx_: Wrapper functions for Scalapack diagonalization routines.
 * - pxxxgvx_post_processing: Post-processing function to handle errors and adjust parameters based on Scalapack's info
 * output.
 * - get_lwork: Helper functions to retrieve the optimal workspace size.
 * - pxxxgvx_diag: Template function to perform diagonalization for different data types.
 *
 * Template instantiations for pxxxgvx_diag are provided for double, float, std::complex<double>, and
 * std::complex<float>.
 */

/**
 * @brief Wrapper function for Scalapack's generalized eigensolver routines.
 * 
 * @param itype Specifies the problem type to be solved.
 * @param jobz Specifies whether to compute eigenvectors.
 * @param range Specifies the range of eigenvalues to be found.
 * @param uplo Specifies whether the upper or lower triangular part of the matrices is referenced.
 * @param n The order of the matrices A and B.
 * @param A The array containing the matrix A.
 * @param ia The row index in the global array A.
 * @param ja The column index in the global array A.
 * @param desca The array descriptor for the distributed matrix A.
 * @param B The array containing the matrix B.
 * @param ib The row index in the global array B.
 * @param jb The column index in the global array B.
 * @param descb The array descriptor for the distributed matrix B.
 * @param vl Lower bound of the interval to be searched for eigenvalues.
 * @param vu Upper bound of the interval to be searched for eigenvalues.
 * @param il Index of the smallest eigenvalue to be returned.
 * @param iu Index of the largest eigenvalue to be returned.
 * @param abstol The absolute error tolerance for the eigenvalues.
 * @param m The total number of eigenvalues found.
 * @param nz The total number of eigenvalues found in the interval (vl, vu].
 * @param w The array to store the eigenvalues.
 * @param orfac The orthogonality factor.
 * @param Z The array to store the eigenvectors.
 * @param iz The row index in the global array Z.
 * @param jz The column index in the global array Z.
 * @param descz The array descriptor for the distributed matrix Z.
 * @param work Workspace array.
 * @param lwork The dimension of the array work.
 * @param rwork Workspace array (not used in this function).
 * @param lrwork The dimension of the array rwork (not used in this function).
 * @param iwork Workspace array.
 * @param liwork The dimension of the array iwork.
 * @param ifail The array to store the indices of the eigenvectors that failed to converge.
 * @param iclustr The array to store the indices of the eigenvalue clusters.
 * @param gap The array to store the gaps between eigenvalue clusters.
 * @param info Output status of the computation.
 * 
 * @note for a uniform interface, rwork and lrwork are input arguments, but not used in pdsygvx_/pssygvx_
 */
void pxxxgvx_(const int* itype,
              const char* jobz,
              const char* range,
              const char* uplo,
              const int* n,
              double* A,
              const int* ia,
              const int* ja,
              const int* desca,
              double* B,
              const int* ib,
              const int* jb,
              const int* descb,
              const double* vl,
              const double* vu,
              const int* il,
              const int* iu,
              const double* abstol,
              int* m,
              int* nz,
              double* w,
              const double* orfac,
              double* Z,
              const int* iz,
              const int* jz,
              const int* descz,
              double* work,
              int* lwork,
              double* rwork, 
              int* lrwork,
              int* iwork,
              int* liwork,
              int* ifail,
              int* iclustr,
              double* gap,
              int* info)
{
    // double
    pdsygvx_(itype,
             jobz,
             range,
             uplo,
             n,
             A,
             ia,
             ja,
             desca,
             B,
             ib,
             jb,
             descb,
             vl,
             vu,
             il,
             iu,
             abstol,
             m,
             nz,
             w,
             orfac,
             Z,
             iz,
             jz,
             descz,
             work,
             lwork,
             iwork,
             liwork,
             ifail,
             iclustr,
             gap,
             info);
}

void pxxxgvx_(const int* itype,
              const char* jobz,
              const char* range,
              const char* uplo,
              const int* n,
              std::complex<double>* A,
              const int* ia,
              const int* ja,
              const int* desca,
              std::complex<double>* B,
              const int* ib,
              const int* jb,
              const int* descb,
              const double* vl,
              const double* vu,
              const int* il,
              const int* iu,
              const double* abstol,
              int* m,
              int* nz,
              double* w,
              const double* orfac,
              std::complex<double>* Z,
              const int* iz,
              const int* jz,
              const int* descz,
              std::complex<double>* work,
              int* lwork,
              double* rwork,
              int* lrwork,
              int* iwork,
              int* liwork,
              int* ifail,
              int* iclustr,
              double* gap,
              int* info)
{
    // std::complex<double>
    pzhegvx_(itype,
             jobz,
             range,
             uplo,
             n,
             A,
             ia,
             ja,
             desca,
             B,
             ib,
             jb,
             descb,
             vl,
             vu,
             il,
             iu,
             abstol,
             m,
             nz,
             w,
             orfac,
             Z,
             iz,
             jz,
             descz,
             work,
             lwork,
             rwork,
             lrwork,
             iwork,
             liwork,
             ifail,
             iclustr,
             gap,
             info);
}

void pxxxgvx_(const int* itype,
              const char* jobz,
              const char* range,
              const char* uplo,
              const int* n,
              float* A,
              const int* ia,
              const int* ja,
              const int* desca,
              float* B,
              const int* ib,
              const int* jb,
              const int* descb,
              const float* vl,
              const float* vu,
              const int* il,
              const int* iu,
              const float* abstol,
              int* m,
              int* nz,
              float* w,
              const float* orfac,
              float* Z,
              const int* iz,
              const int* jz,
              const int* descz,
              float* work,
              int* lwork,
              float* rwork,
              int* lrwork,
              int* iwork,
              int* liwork,
              int* ifail,
              int* iclustr,
              float* gap,
              int* info)
{
    // float
    pssygvx_(itype,
             jobz,
             range,
             uplo,
             n,
             A,
             ia,
             ja,
             desca,
             B,
             ib,
             jb,
             descb,
             vl,
             vu,
             il,
             iu,
             abstol,
             m,
             nz,
             w,
             orfac,
             Z,
             iz,
             jz,
             descz,
             work,
             lwork,
             iwork,
             liwork,
             ifail,
             iclustr,
             gap,
             info);
}

void pxxxgvx_(const int* itype,
              const char* jobz,
              const char* range,
              const char* uplo,
              const int* n,
              std::complex<float>* A,
              const int* ia,
              const int* ja,
              const int* desca,
              std::complex<float>* B,
              const int* ib,
              const int* jb,
              const int* descb,
              const float* vl,
              const float* vu,
              const int* il,
              const int* iu,
              const float* abstol,
              int* m,
              int* nz,
              float* w,
              const float* orfac,
              std::complex<float>* Z,
              const int* iz,
              const int* jz,
              const int* descz,
              std::complex<float>* work,
              int* lwork,
              float* rwork,
              int* lrwork,
              int* iwork,
              int* liwork,
              int* ifail,
              int* iclustr,
              float* gap,
              int* info)
{
    // std::complex<float>
    pchegvx_(itype,
             jobz,
             range,
             uplo,
             n,
             A,
             ia,
             ja,
             desca,
             B,
             ib,
             jb,
             descb,
             vl,
             vu,
             il,
             iu,
             abstol,
             m,
             nz,
             w,
             orfac,
             Z,
             iz,
             jz,
             descz,
             work,
             lwork,
             rwork,
             lrwork,
             iwork,
             liwork,
             ifail,
             iclustr,
             gap,
             info);
}

void pxxxgvx_post_processing(const int info,
                             const std::vector<int>& ifail,
                             const std::vector<int>& iclustr,
                             const int M,
                             const int NZ,
                             const int nbands,
                             int& degeneracy_max)
{
    const std::string str_info = "Scalapack diagonalization: \n    info = " + std::to_string(info) + ".\n";

    if (info == 0)
    {
        return;
    }
    else if (info < 0)
    {
        const int info_negative = -info;
        const std::string str_index = (info_negative > 100)
                                          ? std::to_string(info_negative / 100) + "-th argument "
                                                + std::to_string(info_negative % 100) + "-entry is illegal.\n"
                                          : std::to_string(info_negative) + "-th argument is illegal.\n";
        throw std::runtime_error(str_info + str_index);
    }
    else if (info % 2)
    {
        std::string str_ifail = "ifail = ";
        for (const int i: ifail)
        {
            str_ifail += std::to_string(i) + " ";
        }
        throw std::runtime_error(str_info + str_ifail);
    }
    else if (info / 2 % 2)
    {
        int degeneracy_need = 0;
        for (int irank = 0; irank < iclustr.size() / 2; ++irank)
        {
            degeneracy_need = std::max(degeneracy_need, iclustr[2 * irank + 1] - iclustr[2 * irank]);
        }
        const std::string str_need = "degeneracy_need = " + std::to_string(degeneracy_need) + ".\n";
        const std::string str_saved = "degeneracy_saved = " + std::to_string(degeneracy_max) + ".\n";
        if (degeneracy_need <= degeneracy_max)
        {
            throw std::runtime_error(str_info + str_need + str_saved);
        }
        else
        {
            std::cout << str_need << str_saved;
            degeneracy_max = degeneracy_need;
            return;
        }
    }
    else if (info / 4 % 2)
    {
        const std::string str_M = "M = " + std::to_string(M) + ".\n";
        const std::string str_NZ = "NZ = " + std::to_string(NZ) + ".\n";
        const std::string str_NBANDS = "Number of eigenvalues solved = " + std::to_string(nbands) + ".\n";
        throw std::runtime_error(str_info + str_M + str_NZ + str_NBANDS);
    }
    else if (info / 16 % 2)
    {
        const std::string str_npos = "Not positive definite = " + std::to_string(ifail[0]) + ".\n";
        throw std::runtime_error(str_info + str_npos);
    }
    else
    {
        throw std::runtime_error(str_info);
    }
}

void get_lwork(int& lwork, std::vector<double>& work)
{
    lwork = work[0];
}

void get_lwork(int& lwork, std::vector<float>& work)
{
    lwork = work[0];
}

void get_lwork(int& lwork, std::vector<std::complex<double>>& work)
{
    lwork = work[0].real();
}

void get_lwork(int& lwork, std::vector<std::complex<float>>& work)
{
    lwork = work[0].real();
}

template <typename T>
void pxxxgvx_diag(const int* const desc,
                  const int ncol,
                  const int nrow,
                  const int nbands,
                  const T* const h_mat,
                  const T* const s_mat,
                  typename GetTypeReal<T>::type* const ekb,
                  T* const wfc_2d)
{
    int nprow = 1;
    int npcol = 1;
    int myprow = 0;
    int mypcol = 0;
    Cblacs_gridinfo(desc[1], &nprow, &npcol, &myprow, &mypcol);
    int dsize = nprow * npcol;

    int degeneracy_max = 12; // only used for complex<float> and complex<double>
    while (true)
    {
        std::vector<T> h_tmp(ncol * nrow, 0);
        std::vector<T> s_tmp(ncol * nrow, 0);
        memcpy(h_tmp.data(), h_mat, sizeof(T) * ncol * nrow);
        memcpy(s_tmp.data(), s_mat, sizeof(T) * ncol * nrow);

        int ndim_global = desc[2];
        const char jobz = 'V';
        const char range = 'I';
        const char uplo = 'U';
        const int itype = 1;
        const int il = 1;
        const int iu = nbands;
        const int one = 1;
        int M = 0;
        int NZ = 0;
        int lwork = -1;
        int lrwork = -1;
        int liwork = -1;
        int info = 0;
        const typename GetTypeReal<T>::type abstol = 0;
        const typename GetTypeReal<T>::type orfac = -1;
        const typename GetTypeReal<T>::type vl = 0;
        const typename GetTypeReal<T>::type vu = 0;
        std::vector<T> work(1, 0);
        std::vector<typename GetTypeReal<T>::type> rwork(3, 0); // only used for complex<float> and complex<double>
        std::vector<int> iwork(1, 0);
        std::vector<int> ifail(ndim_global, 0);
        std::vector<int> iclustr(2 * dsize);
        std::vector<typename GetTypeReal<T>::type> gap(dsize);

        pxxxgvx_(&itype,
                 &jobz,
                 &range,
                 &uplo,
                 &ndim_global,
                 h_tmp.data(),
                 &one,
                 &one,
                 desc,
                 s_tmp.data(),
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
                 &lwork,       // is not used for real data type
                 rwork.data(), // is not used for real data type
                 &lrwork,
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);

        if (info)
        {
            throw std::runtime_error("Scalapack diagonalization: \n    info = " + std::to_string(info) + ".\n");
        }

        if (std::is_same<T, std::complex<float>>::value || std::is_same<T, std::complex<double>>::value)
        {
            get_lwork(lwork, work);
            work.resize(lwork, 0);
            liwork = iwork[0];
            iwork.resize(liwork, 0);
            lrwork = rwork[0] + degeneracy_max * ndim_global;
            int maxlrwork = std::max(lrwork, 3);
            rwork.resize(maxlrwork, 0);
        }
        else
        {
            get_lwork(lwork, work);
            work.resize(std::max(lwork, 3), 0);
            liwork = iwork[0];
            iwork.resize(liwork, 0);
        }

        pxxxgvx_(&itype,
                 &jobz,
                 &range,
                 &uplo,
                 &ndim_global,
                 h_tmp.data(),
                 &one,
                 &one,
                 desc,
                 s_tmp.data(),
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
                 rwork.data(), // is not used for real data type
                 &lrwork,      // is not used for real data type
                 iwork.data(),
                 &liwork,
                 ifail.data(),
                 iclustr.data(),
                 gap.data(),
                 &info);

        if (info == 0)
        {
            return;
        }
        pxxxgvx_post_processing(info, ifail, iclustr, M, NZ, nbands, degeneracy_max);

        // break the loop for real data type
        if (std::is_same<T, float>::value || std::is_same<T, double>::value)
        {
            return;
        }
    }
}

// template instantiation
// double
template void pxxxgvx_diag(const int* const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const double* const h_mat,
                           const double* const s_mat,
                           double* const ekb,
                           double* const wfc_2d);

// std::complex<double>                           
template void pxxxgvx_diag(const int* const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const std::complex<double>* const h_mat,
                           const std::complex<double>* const s_mat,
                           double* const ekb,
                           std::complex<double>* const wfc_2d);

// float
template void pxxxgvx_diag(const int* const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const float* const h_mat,
                           const float* const s_mat,
                           float* const ekb,
                           float* const wfc_2d);

// std::complex<float>
template void pxxxgvx_diag(const int* const desc,
                           const int ncol,
                           const int nrow,
                           const int nbands,
                           const std::complex<float>* const h_mat,
                           const std::complex<float>* const s_mat,
                           float* const ekb,
                           std::complex<float>* const wfc_2d);

#endif

} // namespace hsolver