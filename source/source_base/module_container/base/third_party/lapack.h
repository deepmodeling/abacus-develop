/**
 * @file lapack.h
 * @brief This is a direct wrapper of some LAPACK routines.
 * \b Column-Major version.
 * Direct wrapping of standard LAPACK routines. (Column-Major, fortran style)
 *
 * @warning For Row-major version, please refer to \c source/source_base/module_external/lapack_connector.h.
 *
 * @note
 * Some slight modification are made to fit the C++ style for overloading purpose.
 * You can find some function with different parameter list than the original LAPACK routine.
 * And some of these parameters are not referred in the function body. They are included just to
 * ensure the same parameter list for overloaded functions with a uniform name.
 */

#ifndef BASE_THIRD_PARTY_LAPACK_H_
#define BASE_THIRD_PARTY_LAPACK_H_

#include <complex>


#if defined(__CUDA)
#include <base/third_party/cusolver.h>
#elif defined(__ROCM)
#include <base/third_party/hipsolver.h>
#endif

/// This is a wrapper of some LAPACK routines.
/// Direct wrapping of standard LAPACK routines. (column major, fortran style)
/// with some slight modification to fit the C++ style for overloading purpose.

//Naming convention of lapack subroutines : ammxxx, where
//"a" specifies the data type:
//  - d stands for double
//  - z stands for complex double
//"mm" specifies the type of matrix, for example:
//  - he stands for hermitian
//  - sy stands for symmetric
//"xxx" specifies the type of problem, for example:
//  - gv stands for generalized eigenvalue

extern "C"
{
// ILAENV - environment inquiry
int ilaenv_(const int* ispec, const char* name, const char* opts,
            const int* n1, const int* n2, const int* n3, const int* n4);

// Generalized symmetric-definite eigenproblems (divide-and-conquer)
void ssygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             float* a, const int* lda,
             float* b, const int* ldb,
             float* w,
             float* work, const int* lwork,
             int* iwork, const int* liwork,
             int* info);

void dsygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             double* a, const int* lda,
             double* b, const int* ldb,
             double* w,
             double* work, const int* lwork,
             int* iwork, const int* liwork,
             int* info);

void chegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda,
             std::complex<float>* b, const int* ldb,
             float* w,
             std::complex<float>* work, const int* lwork,
             float* rwork, const int* lrwork,
             int* iwork, const int* liwork,
             int* info);

void zhegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<double>* a, const int* lda,
             std::complex<double>* b, const int* ldb,
             double* w,
             std::complex<double>* work, const int* lwork,
             double* rwork, const int* lrwork,
             int* iwork, const int* liwork,
             int* info);

// Generalized symmetric-definite eigenproblems (selected eigenvalues/vectors)
void ssygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
             const int* n, float* A, const int* lda, float* B, const int* ldb,
             const float* vl, const float* vu, const int* il, const int* iu,
             const float* abstol, int* m, float* w, float* Z, const int* ldz,
             float* work, const int* lwork, int* iwork, int* ifail, int* info);

void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
             const int* n, double* A, const int* lda, double* B, const int* ldb,
             const double* vl, const double* vu, const int* il, const int* iu,
             const double* abstol, int* m, double* w, double* Z, const int* ldz,
             double* work, const int* lwork, int* iwork, int* ifail, int* info);

void chegvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
             const int* n, std::complex<float>* A, const int* lda, std::complex<float>* B, const int* ldb,
             const float* vl, const float* vu, const int* il, const int* iu,
             const float* abstol, int* m, float* w, std::complex<float>* Z, const int* ldz,
             std::complex<float>* work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);

void zhegvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
             const int* n, std::complex<double>* A, const int* lda, std::complex<double>* B, const int* ldb,
             const double* vl, const double* vu, const int* il, const int* iu,
             const double* abstol, int* m, double* w, std::complex<double>* Z, const int* ldz,
             std::complex<double>* work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);

// Standard symmetric eigenproblems (selected)
void ssyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             float* a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu,
             const float* abstol, int* m, float* w, float* z, const int* ldz,
             float* work, const int* lwork,  int* iwork, int* ifail, int* info);

void dsyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             double* a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu,
             const double* abstol, int* m, double* w, double* z, const int* ldz,
             double* work, const int* lwork, int* iwork, int* ifail, int* info);

void cheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu,
             const float* abstol, int* m, float* w, std::complex<float>* z, const int* ldz,
             std::complex<float>* work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);

void zheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<double>* a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu,
             const double* abstol, int* m, double* w, std::complex<double>* z, const int* ldz,
             std::complex<double>* work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);

// Standard symmetric eigenproblems (divide-and-conquer)
void ssyevd_(const char* jobz, const char* uplo, const int* n,
             float* a, const int* lda, float* w,
             float* work, const int* lwork,
             int* iwork, const int* liwork, int* info);

void dsyevd_(const char* jobz, const char* uplo, const int* n,
             double* a, const int* lda, double* w,
             double* work, const int* lwork,
             int* iwork, const int* liwork, int* info);

void cheevd_(const char* jobz, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda, float* w,
             std::complex<float>* work, const int* lwork, float* rwork, const int* lrwork,
             int* iwork, const int* liwork, int* info);

void zheevd_(const char* jobz, const char* uplo, const int* n,
             std::complex<double>* a, const int* lda, double* w,
             std::complex<double>* work, const int* lwork, double* rwork, const int* lrwork,
             int* iwork, const int* liwork, int* info);

// Cholesky factorization
void spotrf_(const char* uplo, const int* n, float* A, const int* lda, int* info);
void dpotrf_(const char* uplo, const int* n, double* A, const int* lda, int* info);
void cpotrf_(const char* uplo, const int* n, std::complex<float>* A, const int* lda, int* info);
void zpotrf_(const char* uplo, const int* n, std::complex<double>* A, const int* lda, int* info);

// Inverse using Cholesky factorization
void spotri_(const char* uplo, const int* n, float* A, const int* lda, int* info);
void dpotri_(const char* uplo, const int* n, double* A, const int* lda, int* info);
void cpotri_(const char* uplo, const int* n, std::complex<float>* A, const int* lda, int* info);
void zpotri_(const char* uplo, const int* n, std::complex<double>* A, const int* lda, int* info);

// Inverse of triangular matrix
void strtri_(const char* uplo, const char* diag, const int* n, float* a, const int* lda, int* info);
void dtrtri_(const char* uplo, const char* diag, const int* n, double* a, const int* lda, int* info);
void ctrtri_(const char* uplo, const char* diag, const int* n, std::complex<float>* a, const int* lda, int* info);
void ztrtri_(const char* uplo, const char* diag, const int* n, std::complex<double>* a, const int* lda, int* info);

// LU factorization
void sgetrf_(const int* m, const int* n, float* a, const int* lda, int* ipiv, int* info);
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info);
void cgetrf_(const int* m, const int* n, std::complex<float>* a, const int* lda, int* ipiv, int* info);
void zgetrf_(const int* m, const int* n, std::complex<double>* a, const int* lda, int* ipiv, int* info);

// Inverse using LU factorization
void sgetri_(const int* n, float* A, const int* lda, const int* ipiv, float* work, const int* lwork, int* info);
void dgetri_(const int* n, double* A, const int* lda, const int* ipiv, double* work, const int* lwork, int* info);
void cgetri_(const int* n, std::complex<float>* A, const int* lda, const int* ipiv, std::complex<float>* work, const int* lwork, int* info);
void zgetri_(const int* n, std::complex<double>* A, const int* lda, const int* ipiv, std::complex<double>* work, const int* lwork, int* info);

// Solve linear system using LU factorization
void sgetrs_(const char* trans, const int* n, const int* nrhs,
             const float* A, const int* lda, const int* ipiv,
             float* B, const int* ldb, int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs,
             const double* A, const int* lda, const int* ipiv,
             double* B, const int* ldb, int* info);
void cgetrs_(const char* trans, const int* n, const int* nrhs,
             const std::complex<float>* A, const int* lda, const int* ipiv,
             std::complex<float>* B, const int* ldb, int* info);
void zgetrs_(const char* trans, const int* n, const int* nrhs,
             const std::complex<double>* A, const int* lda, const int* ipiv,
             std::complex<double>* B, const int* ldb, int* info);
}

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
namespace container {
namespace lapackConnector
{
static inline
int ilaenv( int ispec, const char *name,const char *opts,const int n1,const int n2,
            const int n3,const int n4)
{
    const int nb = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
    return nb;
}
// wrap function of fortran lapack routine zhegvd. (pointer version)
static inline
void hegvd(const int itype, const char jobz, const char uplo, const int n,
            float* a, const int lda,
            float* b, const int ldb, float* w,
            float* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    ssygvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork,
            iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void hegvd(const int itype, const char jobz, const char uplo, const int n,
            double* a, const int lda,
            double* b, const int ldb, double* w,
            double* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    dsygvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork,
            iwork, &liwork, &info);
}
static inline
void hegvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<float>* a, const int lda,
            std::complex<float>* b, const int ldb, float* w,
            std::complex<float>* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    chegvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void hegvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<double>* a, const int lda,
            std::complex<double>* b, const int ldb, double* w,
            std::complex<double>* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    zhegvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
}

// Note
// rwork is only needed for complex version
// and we include rwork in the function parameter list
// for simplicity of function overloading
// and unification of function parameter list
static inline
void hegvx(const int itype, const char jobz, const char range, const char uplo, const int n,
            float* a, const int lda, float* b, const int ldb,
            const float vl, const float vu, const int il, const int iu, const float abstol,
            int m, float* w, float* z, const int ldz,
            float* work, const int lwork, float* rwork, int* iwork, int* ifail, int& info)
{
    ssygvx_(&itype, &jobz, &range, &uplo, &n,
            a, &lda, b, &ldb,
            &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, iwork, ifail, &info);
}

static inline
void hegvx(const int itype, const char jobz, const char range, const char uplo, const int n,
            double* a, const int lda, double* b, const int ldb,
            const double vl, const double vu, const int il, const int iu, const double abstol,
            int m, double* w, double* z, const int ldz,
            double* work, const int lwork, double* rwork, int* iwork, int* ifail, int& info)
{
    dsygvx_(&itype, &jobz, &range, &uplo, &n,
            a, &lda, b, &ldb,
            &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, iwork, ifail, &info);
}

static inline
void hegvx(const int itype, const char jobz, const char range, const char uplo, const int n,
            std::complex<float>* a, const int lda, std::complex<float>* b, const int ldb,
            const float vl, const float vu, const int il, const int iu, const float abstol,
            int m, float* w, std::complex<float>* z, const int ldz,
            std::complex<float>* work, const int lwork, float* rwork, int* iwork, int* ifail, int& info)
{
    chegvx_(&itype, &jobz, &range, &uplo, &n,
            a, &lda, b, &ldb,
            &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}

static inline
void hegvx(const int itype, const char jobz, const char range, const char uplo, const int n,
            std::complex<double>* a, const int lda, std::complex<double>* b, const int ldb,
            const double vl, const double vu, const int il, const int iu, const double abstol,
            int m, double* w, std::complex<double>* z, const int ldz,
            std::complex<double>* work, const int lwork, double* rwork, int* iwork, int* ifail, int& info)
{
    zhegvx_(&itype, &jobz, &range, &uplo, &n,
            a, &lda, b, &ldb,
            &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}


// wrap function of fortran lapack routine zheevx.
static inline
void heevx(const char jobz, const char range, const char uplo, const int n,
             float* a, const int lda,
             const float vl, const float vu, const int il, const int iu, const float abstol,
             int m, float* w, float* z, const int ldz,
             float* work, const int lwork, float* rwork, int* iwork, int* ifail, int info)
{
    ssyevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork,  iwork, ifail, &info);
}
// wrap function of fortran lapack routine zheevx.
static inline
void heevx(const char jobz, const char range, const char uplo, const int n,
            double* a, const int lda,
            const double vl, const double vu, const int il, const int iu, const double abstol,
            int m, double* w, double* z, const int ldz,
            double* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
{
    dsyevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, iwork, ifail, &info);
}
static inline
void heevx(const char jobz, const char range, const char uplo, const int n,
             std::complex<float>* a, const int lda,
             const float vl, const float vu, const int il, const int iu, const float abstol,
             int m, float* w, std::complex<float>* z, const int ldz,
             std::complex<float>* work, const int lwork, float* rwork, int* iwork, int* ifail, int info)
{
    cheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}
// wrap function of fortran lapack routine zheevx.
static inline
void heevx(const char jobz, const char range, const char uplo, const int n,
             std::complex<double>* a, const int lda,
             const double vl, const double vu, const int il, const int iu, const double abstol,
             int m, double* w, std::complex<double>* z, const int ldz,
             std::complex<double>* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
{
    zheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}

static inline
void heevd(const char jobz, const char uplo, const int n,
            float* a, const int lda, float* w,
            float* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    ssyevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork,
             iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void heevd(const char jobz, const char uplo, const int n,
            double* a, const int lda, double* w,
            double* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    dsyevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork,
             iwork, &liwork, &info);
}
static inline
void heevd(const char jobz, const char uplo, const int n,
            std::complex<float>* a, const int lda, float* w,
            std::complex<float>* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    cheevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork, rwork, &lrwork,
             iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void heevd(const char jobz, const char uplo, const int n,
            std::complex<double>* a, const int lda, double* w,
            std::complex<double>* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    zheevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork, rwork, &lrwork,
             iwork, &liwork, &info);
}

static inline
void potrf( const char &uplo, const int &n, float* A, const int &lda, int &info )
{
	spotrf_(&uplo, &n, A, &lda, &info );
}
static inline
void potrf( const char &uplo, const int &n, double* A, const int &lda, int &info )
{
	dpotrf_(&uplo, &n, A, &lda, &info );
}
static inline
void potrf( const char &uplo, const int &n, std::complex<float>* A, const int &lda, int &info )
{
	cpotrf_(&uplo, &n, A, &lda, &info );
}
static inline
void potrf( const char &uplo, const int &n, std::complex<double>* A, const int &lda, int &info )
{
	zpotrf_( &uplo, &n, A, &lda, &info );
}

static inline
void trtri( const char &uplo, const char &diag, const int &n, float* A, const int &lda, int &info )
{
    strtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, double* A, const int &lda, int &info)
{
    dtrtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, std::complex<float>* A, const int &lda, int &info )
{
    ctrtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, std::complex<double>* A, const int &lda, int &info)
{
    ztrtri_( &uplo, &diag, &n, A, &lda, &info);
}

static inline
void getrf(const int m, const int n, float* A, const int lda, int* ipiv, int &info)
{
    sgetrf_(&m, &n, A, &lda, ipiv, &info);
}
static inline
void getrf(const int m, const int n, double* A, const int lda, int* ipiv, int &info)
{
    dgetrf_(&m, &n, A, &lda, ipiv, &info);
}
static inline
void getrf(const int m, const int n, std::complex<float>* A, const int lda, int* ipiv, int &info)
{
    cgetrf_(&m, &n, A, &lda, ipiv, &info);
}
static inline
void getrf(const int m, const int n, std::complex<double>* A, const int lda, int* ipiv, int &info)
{
    zgetrf_(&m, &n, A, &lda, ipiv, &info);
}

static inline
void getri(const int n, float* A, const int lda, const int* ipiv, float* work, const int lwork, int& info)
{
    sgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
}
static inline
void getri(const int n, double* A, const int lda, const int* ipiv, double* work, const int lwork, int& info)
{
    dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
}
static inline
void getri(const int n, std::complex<float>* A, const int lda, const int* ipiv, std::complex<float>* work, const int lwork, int& info)
{
    cgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
}
static inline
void getri(const int n, std::complex<double>* A, const int lda, const int* ipiv, std::complex<double>* work, const int lwork, int& info)
{
    zgetri_(&n, A, &lda, ipiv, work, &lwork, &info);
}

static inline
void getrs(const char& trans, const int n, const int nrhs, float* A, const int lda, const int* ipiv, float* B, const int ldb, int& info)
{
    sgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
}
static inline
void getrs(const char& trans, const int n, const int nrhs, double* A, const int lda, const int* ipiv, double* B, const int ldb, int& info)
{
    dgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
}
static inline
void getrs(const char& trans, const int n, const int nrhs, std::complex<float>* A, const int lda, const int* ipiv, std::complex<float>* B, const int ldb, int& info)
{
    cgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
}
static inline
void getrs(const char& trans, const int n, const int nrhs, std::complex<double>* A, const int lda, const int* ipiv, std::complex<double>* B, const int ldb, int& info)
{
    zgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
}

} // namespace lapackConnector
} // namespace container

#endif  // BASE_THIRD_PARTY_LAPACK_H_
