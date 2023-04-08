// This file defines the math lib wrapper for output information before executing computations.

/*
When INFO is defined in cmake configure, a macro of GATHER_INFO will be defined.
This macro will be used to output information before executing computations.
Results will output to OUT/math_info.log, see ModuleBase::Global_File::make_dir_out .
*/

#include <ostream>
#undef GATHER_INFO
#include "module_base/blas_connector.h"
#include "module_base/lapack_connector.h"

#include <iostream>

void zgemm_i(const char *transa,
             const char *transb,
             const int *m,
             const int *n,
             const int *k,
             const std::complex<double> *alpha,
             const std::complex<double> *a,
             const int *lda,
             const std::complex<double> *b,
             const int *ldb,
             const std::complex<double> *beta,
             std::complex<double> *c,
             const int *ldc)
{
    GlobalV::ofs_info.unsetf(std::ios_base::floatfield);
    GlobalV::ofs_info << "zgemm " << *transa << " " << *transb << " " << *m << " " << *n << " "
                      << *k << " " << *alpha << " " << *lda << " " << *ldb << " " << *beta << " " << *ldc << std::endl;
    zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

void zaxpy_i(const int *N,
             const std::complex<double> *alpha,
             const std::complex<double> *X,
             const int *incX,
             std::complex<double> *Y,
             const int *incY)
{
    // std::cout << "zaxpy " << *N << std::endl;
    // alpha is a coefficient
    // incX, incY is always 1
    zaxpy_(N, alpha, X, incX, Y, incY);
}

void zhegvx_i(const int *itype,
              const char *jobz,
              const char *range,
              const char *uplo,
              const int *n,
              std::complex<double> *a,
              const int *lda,
              std::complex<double> *b,
              const int *ldb,
              const double *vl,
              const double *vu,
              const int *il,
              const int *iu,
              const double *abstol,
              const int *m,
              double *w,
              std::complex<double> *z,
              const int *ldz,
              std::complex<double> *work,
              const int *lwork,
              double *rwork,
              int *iwork,
              int *ifail,
              int *info)
{
    GlobalV::ofs_info.unsetf(std::ios_base::floatfield);
    GlobalV::ofs_info << "zhegvx " << *itype << " " << *jobz << " " << *range << " " << *uplo
                      << " " << *n << " " << *lda << " " << *ldb << " " << *vl << " " << *vu << " " << *il << " " << *iu
                      << " " << *abstol << " " << *m << " " << *lwork << " " << *info << std::endl;
    zhegvx_(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork,
            iwork, ifail, info);
}
