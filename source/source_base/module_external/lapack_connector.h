#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>
#include <stdexcept>
#include <iostream>
#include <cassert>
#include "../matrix.h"
#include "../complexmatrix.h"
#include "../global_function.h"

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
namespace LapackConnector
{
enum MatrixLayout
{
    RowMajor,
    ColMajor
};

void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<float>* a, int lda, std::complex<float>* b, int ldb, float* w);
void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<double>* a, int lda, std::complex<double>* b, int ldb, double* w);
void hegv(MatrixLayout layout, int itype, char jobz, char uplo, int n, double* a, int lda, double* b, int ldb, double* w);
void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, float* a, int lda, float* b, int ldb, float* w);
void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, double* a, int lda, double* b, int ldb, double* w);
void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<float>* a, int lda, std::complex<float>* b, int ldb, float* w);
void hegvd(MatrixLayout layout, int itype, char jobz, char uplo, int n, std::complex<double>* a, int lda, std::complex<double>* b, int ldb, double* w);
void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           std::complex<float>* a, int lda, std::complex<float>* b, int ldb,
           float vl, float vu, int il, int iu, float abstol, int* m,
           float* w, std::complex<float>* z, int ldz, int* ifail);
void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           std::complex<double>* a, int lda, std::complex<double>* b, int ldb,
           double vl, double vu, int il, int iu, double abstol, int* m,
           double* w, std::complex<double>* z, int ldz, int* ifail);
void hegvx(MatrixLayout layout, int itype, char jobz, char range, char uplo, int n,
           double* a, int lda, double* b, int ldb,
           double vl, double vu, int il, int iu, double abstol, int* m,
           double* w, double* z, int ldz, int* ifail);
void getrf(MatrixLayout layout, int m, int n, std::complex<double>* a, int lda, int* ipiv);
void getri(MatrixLayout layout, std::complex<double>* a, int lda, const int* ipiv);

void potrf(MatrixLayout layout, char uplo, int n, float* a, int lda);
void potrf(MatrixLayout layout, char uplo, int n, double* a, int lda);
void potrf(MatrixLayout layout, char uplo, int n, std::complex<float>* a, int lda);
void potrf(MatrixLayout layout, char uplo, int n, std::complex<double>* a, int lda);

void potri(MatrixLayout layout, char uplo, int n, float* a, int lda);
void potri(MatrixLayout layout, char uplo, int n, double* a, int lda);
void potri(MatrixLayout layout, char uplo, int n, std::complex<float>* a, int lda);
void potri(MatrixLayout layout, char uplo, int n, std::complex<double>* a, int lda);

void heev(MatrixLayout layout, char jobz, char uplo, int n, std::complex<float>* a, int lda, float* w);
void heev(MatrixLayout layout, char jobz, char uplo, int n, std::complex<double>* a, int lda, double* w);
void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, float* a, int lda, float vl,
           float vu, int il, int iu, float abstol, int* m, float* w, float* z, int ldz, int* ifail);
void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, double* a, int lda, double vl,
           double vu, int il, int iu, double abstol, int* m, double* w, double* z, int ldz, int* ifail);
void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, std::complex<float>* a, int lda,
           float vl, float vu, int il, int iu, float abstol, int* m, float* w,
           std::complex<float>* z, int ldz, int* ifail);
void heevx(MatrixLayout layout, char jobz, char range, char uplo, int n, std::complex<double>* a, int lda,
           double vl, double vu, int il, int iu, double abstol, int* m, double* w,
           std::complex<double>* z, int ldz, int* ifail);
void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           float* a, int lda, float* w);
void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           double* a, int lda, double* w);
void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           std::complex<float>* a, int lda, float* w);
void heevd(MatrixLayout layout, char jobz, char uplo, int n,
           std::complex<double>* a, int lda, double* w);
void syev(MatrixLayout layout, char jobz, char uplo, int n, double* a, int lda, double* w);

void geev(MatrixLayout layout, char jobvl, char jobvr, int n, double* a, int lda,
          double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr);
void geev(MatrixLayout layout, char jobvl, char jobvr, int n, std::complex<double>* a, int lda,
          std::complex<double>* w, std::complex<double>* vl, int ldvl, std::complex<double>* vr, int ldvr);

void getrf(MatrixLayout layout, int m, int n, float* a, int lda, int* ipiv);
void getrf(MatrixLayout layout, int m, int n, double* a, int lda, int* ipiv);
void getrf(MatrixLayout layout, int m, int n, std::complex<float>* a, int lda, int* ipiv);
void getrf(MatrixLayout layout, int m, int n, std::complex<double>* a, int lda, int* ipiv);
void getri(MatrixLayout layout, int n, float* a, int lda, const int* ipiv);
void getri(MatrixLayout layout, int n, double* a, int lda, const int* ipiv);
void getri(MatrixLayout layout, int n, std::complex<float>* a, int lda, const int* ipiv);
void getri(MatrixLayout layout, int n, std::complex<double>* a, int lda, const int* ipiv);
void getrs(MatrixLayout layout, char trans, int n, int nrhs, const float* a, int lda, const int* ipiv, float* b, int ldb);
void getrs(MatrixLayout layout, char trans, int n, int nrhs, const double* a, int lda, const int* ipiv, double* b, int ldb);
void getrs(MatrixLayout layout, char trans, int n, int nrhs, const std::complex<float>* a, int lda, const int* ipiv, std::complex<float>* b, int ldb);
void getrs(MatrixLayout layout, char trans, int n, int nrhs, const std::complex<double>* a, int lda, const int* ipiv, std::complex<double>* b, int ldb);
void sytrf(MatrixLayout layout, char uplo, int n, double* a, int lda, int* ipiv);
void sytri(MatrixLayout layout, char uplo, int n, double* a, int lda, const int* ipiv);

void gtsv(MatrixLayout layout, int n, int nrhs, double* dl, double* d, double* du, double* b, int ldb);
void sysv(MatrixLayout layout, char uplo, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);

void trtri(MatrixLayout layout, char uplo, char diag, int n, float* a, int lda);
void trtri(MatrixLayout layout, char uplo, char diag, int n, double* a, int lda);
void trtri(MatrixLayout layout, char uplo, char diag, int n, std::complex<double>* a, int lda);
void trtri(MatrixLayout layout, char uplo, char diag, int n, std::complex<float>* a, int lda);
}
#endif  // LAPACKCONNECTOR_HPP
