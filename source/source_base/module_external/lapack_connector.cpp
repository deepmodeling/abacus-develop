#include <lapacke.h>
#include "lapack_connector.h"

namespace LapackConnector
{
int toLapackLayout(MatrixLayout layout)
{
    return (layout == MatrixLayout::RowMajor) ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
}

void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<float>* a, int lda, std::complex<float>* b, int ldb, float* w)
{
    
    int info = LAPACKE_chegv(toLapackLayout(layout), itype, jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, reinterpret_cast<lapack_complex_float*>(b), ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK chegv failed with info = " + std::to_string(info));
    }
}

void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<double>* a, int lda, std::complex<double>* b, int ldb, double* w)
{
    int info = LAPACKE_zhegv(toLapackLayout(layout), itype, jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, reinterpret_cast<lapack_complex_double*>(b), ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zhegv failed with info = " + std::to_string(info));
    }
}

void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, double* a, int lda, double* b, int ldb, double* w)
{
    int info = LAPACKE_dsygv(toLapackLayout(layout), itype, jobz, uplo, n, a, lda, b, ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsygv failed with info = " + std::to_string(info));
    }
}

void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, float* a, int lda, float* b, int ldb, float* w)
{
    int info = LAPACKE_ssygvd(toLapackLayout(layout), itype, jobz, uplo, n, a, lda, b, ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK ssygvd failed with info = " + std::to_string(info));
    }
}

void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, double* a, int lda, double* b, int ldb, double* w)
{
    int info = LAPACKE_dsygvd(toLapackLayout(layout), itype, jobz, uplo, n, a, lda, b, ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsygvd failed with info = " + std::to_string(info));
    }
}

void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<float>* a, int lda, std::complex<float>* b, int ldb, float* w)
{
    int info = LAPACKE_chegvd(toLapackLayout(layout), itype, jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, reinterpret_cast<lapack_complex_float*>(b), ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK chegvd failed with info = " + std::to_string(info));
    }
}

void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<double>* a, int lda, std::complex<double>* b, int ldb, double* w)
{
    int info = LAPACKE_zhegvd(toLapackLayout(layout), itype, jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, reinterpret_cast<lapack_complex_double*>(b), ldb, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zhegvd failed with info = " + std::to_string(info));
    }
}

void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           std::complex<float>* a, int lda, std::complex<float>* b, int ldb,
           float vl, float vu, int il, int iu, float abstol, int* m,
           float* w, std::complex<float>* z, int ldz, int* ifail)
{
    int info = LAPACKE_chegvx(toLapackLayout(layout), itype, jobz, range, uplo, n,
                             reinterpret_cast<lapack_complex_float*>(a), lda,
                             reinterpret_cast<lapack_complex_float*>(b), ldb,
                             vl, vu, il, iu, abstol, m, w,
                             reinterpret_cast<lapack_complex_float*>(z), ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK chegvx failed with info = " + std::to_string(info));
    }
}

void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           std::complex<double>* a, int lda, std::complex<double>* b, int ldb,
           double vl, double vu, int il, int iu, double abstol, int* m,
           double* w, std::complex<double>* z, int ldz, int* ifail)
{
    int info = LAPACKE_zhegvx(toLapackLayout(layout), itype, jobz, range, uplo, n,
                             reinterpret_cast<lapack_complex_double*>(a), lda,
                             reinterpret_cast<lapack_complex_double*>(b), ldb,
                             vl, vu, il, iu, abstol, m, w,
                             reinterpret_cast<lapack_complex_double*>(z), ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zhegvx failed with info = " + std::to_string(info));
    }
}

void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           double* a, int lda, double* b, int ldb,
           double vl, double vu, int il, int iu, double abstol, int* m,
           double* w, double* z, int ldz, int* ifail)
{
    int info = LAPACKE_dsygvx(toLapackLayout(layout), itype, jobz, range, uplo, n,
                             a, lda, b, ldb,
                             vl, vu, il, iu, abstol, m, w,
                             z, ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsygvx failed with info = " + std::to_string(info));
    }
}

void potrf(MatrixLayout layout, char uplo, int n, float* a, int lda)
{
    int info = LAPACKE_spotrf(toLapackLayout(layout), uplo, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK spotrf failed with info = " + std::to_string(info));
    }
}

void potrf(MatrixLayout layout, char uplo, int n, double* a, int lda)
{
    int info = LAPACKE_dpotrf(toLapackLayout(layout), uplo, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dpotrf failed with info = " + std::to_string(info));
    }
}

void potrf(MatrixLayout layout, char uplo, int n, std::complex<float>* a, int lda)
{
    int info = LAPACKE_cpotrf(toLapackLayout(layout), uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cpotrf failed with info = " + std::to_string(info));
    }
}

void potrf(MatrixLayout layout, char uplo, int n, std::complex<double>* a, int lda)
{
    int info = LAPACKE_zpotrf(toLapackLayout(layout), uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zpotrf failed with info = " + std::to_string(info));
    }
}

void potri(MatrixLayout layout, char uplo, int n, float* a, int lda)
{
    int info = LAPACKE_spotri(toLapackLayout(layout), uplo, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK spotri failed with info = " + std::to_string(info));
    }
}

void potri(MatrixLayout layout, char uplo, int n, double* a, int lda)
{
    int info = LAPACKE_dpotri(toLapackLayout(layout), uplo, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dpotri failed with info = " + std::to_string(info));
    }
}

void potri(MatrixLayout layout, char uplo, int n, std::complex<float>* a, int lda)
{
    int info = LAPACKE_cpotri(toLapackLayout(layout), uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cpotri failed with info = " + std::to_string(info));
    }
}

void potri(MatrixLayout layout, char uplo, int n, std::complex<double>* a, int lda)
{
    int info = LAPACKE_zpotri(toLapackLayout(layout), uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zpotri failed with info = " + std::to_string(info));
    }
}

void heev(MatrixLayout layout, char jobz, char uplo, int n, std::complex<float>* a, int lda, float* w)
{
    int info = LAPACKE_cheev(toLapackLayout(layout), jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cheev failed with info = " + std::to_string(info));
    }
}

void heev(MatrixLayout layout, char jobz, char uplo, int n, std::complex<double>* a, int lda, double* w)
{
    int info = LAPACKE_zheev(toLapackLayout(layout), jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zheev failed with info = " + std::to_string(info));
    }
} 

void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, float* a, int lda, float vl,
           float vu, int il, int iu, float abstol, int* m, float* w, float* z, int ldz, int* ifail)
{
    int info = LAPACKE_ssyevx(toLapackLayout(layout), jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK ssyevx failed with info = " + std::to_string(info));
    }
}

void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, double* a, int lda, double vl,
           double vu, int il, int iu, double abstol, int* m, double* w, double* z, int ldz, int* ifail)
{
    int info = LAPACKE_dsyevx(toLapackLayout(layout), jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsyevx failed with info = " + std::to_string(info));
    }
}

void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, std::complex<float>* a, int lda,
           float vl, float vu, int il, int iu, float abstol, int* m, float* w,
           std::complex<float>* z, int ldz, int* ifail)
{
    int info = LAPACKE_cheevx(toLapackLayout(layout), jobz, range, uplo, n,
                             reinterpret_cast<lapack_complex_float*>(a), lda,
                             vl, vu, il, iu, abstol, m, w,
                             reinterpret_cast<lapack_complex_float*>(z), ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cheevx failed with info = " + std::to_string(info));
    }
}

void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, std::complex<double>* a, int lda,
           double vl, double vu, int il, int iu, double abstol, int* m, double* w,
           std::complex<double>* z, int ldz, int* ifail)
{
    int info = LAPACKE_zheevx(toLapackLayout(layout), jobz, range, uplo, n,
                             reinterpret_cast<lapack_complex_double*>(a), lda,
                             vl, vu, il, iu, abstol, m, w,
                             reinterpret_cast<lapack_complex_double*>(z), ldz, ifail);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zheevx failed with info = " + std::to_string(info));
    }
}

void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           float* a, int lda, float* w)
{
    int info = LAPACKE_ssyevd(toLapackLayout(layout), jobz, uplo, n, a, lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK ssyevd failed with info = " + std::to_string(info));
    }
}

void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           double* a, int lda, double* w)
{
    int info = LAPACKE_dsyevd(toLapackLayout(layout), jobz, uplo, n, a, lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsyevd failed with info = " + std::to_string(info));
    }
}

void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           std::complex<float>* a, int lda, float* w)
{
    int info = LAPACKE_cheevd(toLapackLayout(layout), jobz, uplo, n, reinterpret_cast<lapack_complex_float*>(a), lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cheevd failed with info = " + std::to_string(info));
    }
}

void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           std::complex<double>* a, int lda, double* w)
{
    int info = LAPACKE_zheevd(toLapackLayout(layout), jobz, uplo, n, reinterpret_cast<lapack_complex_double*>(a), lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zheevd failed with info = " + std::to_string(info));
    }
}

void syev(MatrixLayout layout, char jobz, char uplo, int n, double* a, int lda, double* w)
{
    int info = LAPACKE_dsyev(toLapackLayout(layout), jobz, uplo, n, a, lda, w);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsyev failed with info = " + std::to_string(info));
    }
}

void geev(MatrixLayout layout, char jobvl, char jobvr, int n, double* a, int lda,
          double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr)
{
    int info = LAPACKE_dgeev(toLapackLayout(layout), jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dgeev failed with info = " + std::to_string(info));
    }
}

void geev(MatrixLayout layout, char jobvl, char jobvr, int n, std::complex<double>* a, int lda,
          std::complex<double>* w, std::complex<double>* vl, int ldvl, std::complex<double>* vr, int ldvr)
{
    int info = LAPACKE_zgeev(toLapackLayout(layout), jobvl, jobvr, n, reinterpret_cast<lapack_complex_double*>(a), lda,
                             reinterpret_cast<lapack_complex_double*>(w), reinterpret_cast<lapack_complex_double*>(vl), ldvl,
                             reinterpret_cast<lapack_complex_double*>(vr), ldvr);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zgeev failed with info = " + std::to_string(info));
    }
}

void getrf(MatrixLayout layout, int m, int n, float* a, int lda, int* ipiv)
{
    int info = LAPACKE_sgetrf(toLapackLayout(layout), m, n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK sgetrf failed with info = " + std::to_string(info));
    }
}

void getrf(MatrixLayout layout, int m, int n, double* a, int lda, int* ipiv)
{
    int info = LAPACKE_dgetrf(toLapackLayout(layout), m, n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dgetrf failed with info = " + std::to_string(info));
    }
}

void getrf(MatrixLayout layout, int m, int n, std::complex<float>* a, int lda, int* ipiv)
{
    int info = LAPACKE_cgetrf(toLapackLayout(layout), m, n, reinterpret_cast<lapack_complex_float*>(a), lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cgetrf failed with info = " + std::to_string(info));
    }
}

void getrf(MatrixLayout layout, int m, int n, std::complex<double>* a, int lda, int* ipiv)
{
    int info = LAPACKE_zgetrf(toLapackLayout(layout), m, n, reinterpret_cast<lapack_complex_double*>(a), lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zgetrf failed with info = " + std::to_string(info));
    }
}

void getri(MatrixLayout layout, int n, float* a, int lda, const int* ipiv)
{
    int info = LAPACKE_sgetri(toLapackLayout(layout), n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK sgetri failed with info = " + std::to_string(info));
    }
}

void getri(MatrixLayout layout, int n, double* a, int lda, const int* ipiv)
{
    int info = LAPACKE_dgetri(toLapackLayout(layout), n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dgetri failed with info = " + std::to_string(info));
    }
}

void getri(MatrixLayout layout, int n, std::complex<float>* a, int lda, const int* ipiv)
{
    int info = LAPACKE_cgetri(toLapackLayout(layout), n, reinterpret_cast<lapack_complex_float*>(a), lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cgetri failed with info = " + std::to_string(info));
    }
}

void getri(MatrixLayout layout, int n, std::complex<double>* a, int lda, const int* ipiv)
{
    int info = LAPACKE_zgetri(toLapackLayout(layout), n, reinterpret_cast<lapack_complex_double*>(a), lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zgetri failed with info = " + std::to_string(info));
    }
}

void getrs(MatrixLayout layout, char trans, int n, int nrhs, const float* a, int lda, const int* ipiv, float* b, int ldb)
{
    int info = LAPACKE_sgetrs(toLapackLayout(layout), trans, n, nrhs, a, lda, ipiv, b, ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK sgetrs failed with info = " + std::to_string(info));
    }
}

void getrs(MatrixLayout layout, char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb)
{
    int info = LAPACKE_dgetrs(toLapackLayout(layout), trans, n, nrhs, a, lda, ipiv, b, ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dgetrs failed with info = " + std::to_string(info));
    }
}

void getrs(MatrixLayout layout, char trans, int n, int nrhs, const std::complex<float>* a, int lda, const int* ipiv, std::complex<float>* b, int ldb)
{
    int info = LAPACKE_cgetrs(toLapackLayout(layout), trans, n, nrhs, reinterpret_cast<const lapack_complex_float*>(a), lda, ipiv, reinterpret_cast<lapack_complex_float*>(b), ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK cgetrs failed with info = " + std::to_string(info));
    }
}

void getrs(MatrixLayout layout, char trans, int n, int nrhs, const std::complex<double>* a, int lda, const int* ipiv, std::complex<double>* b, int ldb)
{
    int info = LAPACKE_zgetrs(toLapackLayout(layout), trans, n, nrhs, reinterpret_cast<const lapack_complex_double*>(a), lda, ipiv, reinterpret_cast<lapack_complex_double*>(b), ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK zgetrs failed with info = " + std::to_string(info));
    }
}

void sytrf(MatrixLayout layout, char uplo, int n, double* a, int lda, int* ipiv)
{
    int info = LAPACKE_dsytrf(toLapackLayout(layout), uplo, n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsytrf failed with info = " + std::to_string(info));
    }
}

void sytri(MatrixLayout layout, char uplo, int n, double* a, int lda, const int* ipiv)
{
    int info = LAPACKE_dsytri(toLapackLayout(layout), uplo, n, a, lda, ipiv);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsytri failed with info = " + std::to_string(info));
    }
}

void trtri(MatrixLayout layout, char uplo, char diag, int n, float* a, int lda)
{
    int info = LAPACKE_strtri(toLapackLayout(layout), uplo, diag, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK strtri failed with info = " + std::to_string(info));
    }
}

void trtri(MatrixLayout layout, char uplo, char diag, int n, double* a, int lda)
{
    int info = LAPACKE_dtrtri(toLapackLayout(layout), uplo, diag, n, a, lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dtrtri failed with info = " + std::to_string(info));
    }
}

void trtri(MatrixLayout layout, char uplo, char diag, int n, std::complex<double>* a, int lda)
{
    int info = LAPACKE_ztrtri(toLapackLayout(layout), uplo, diag, n, reinterpret_cast<lapack_complex_double*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK ztrtri failed with info = " + std::to_string(info));
    }
}

void trtri(MatrixLayout layout, char uplo, char diag, int n, std::complex<float>* a, int lda)
{
    int info = LAPACKE_ctrtri(toLapackLayout(layout), uplo, diag, n, reinterpret_cast<lapack_complex_float*>(a), lda);
    if (info != 0) {
        throw std::invalid_argument("LAPACK ctrtri failed with info = " + std::to_string(info));
    }
}

void gtsv(MatrixLayout layout, int n, int nrhs, double* dl, double* d, double* du, double* b, int ldb)
{
    int info = LAPACKE_dgtsv(toLapackLayout(layout), n, nrhs, dl, d, du, b, ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dgtsv failed with info = " + std::to_string(info));
    }
}

void sysv(MatrixLayout layout, char uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb)
{
    int info = LAPACKE_dsysv(toLapackLayout(layout), uplo, n, nrhs, a, lda, ipiv, b, ldb);
    if (info != 0) {
        throw std::invalid_argument("LAPACK dsysv failed with info = " + std::to_string(info));
    }
}
}  // namespace LapackConnector

