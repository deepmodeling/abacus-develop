#ifndef BASE_THIRD_PARTY_LAPACK_H_
#define BASE_THIRD_PARTY_LAPACK_H_

#include <complex>
#include "source_base/macros.h"
#include "source_base/module_external/lapack_connector.h"

#if defined(__CUDA)
#include <base/third_party/cusolver.h>
#elif defined(__ROCM)
#include <base/third_party/hipsolver.h>
#endif

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
namespace container {
namespace lapackConnector
{

template<typename T>
using Real = typename GetTypeReal<T>::type;
// wrap function of fortran lapack routine zhegvd. (pointer version)
template<typename T>
void dngvd(int itype, char jobz, char uplo, int n,
            T* a, int lda,
            T* b, int ldb, Real<T>* w)
{
    LapackConnector::hegvd(
        LapackConnector::ColMajor, itype, jobz, uplo, n,
        a, lda, b, ldb, w);
}

// wrap function of fortran lapack routine zheevx.
template<typename T>
void dnevx( int itype, char jobz, char range, char uplo, int n,
             T* a, int lda,
             Real<T> vl, Real<T> vu, int il, int iu, Real<T> abstol,
             int m, Real<T>* w, T* z, int ldz, int* ifail)
{
    LapackConnector::heevx(
        LapackConnector::ColMajor, jobz, range, uplo, n,
        a, lda, vl, vu, il, iu,
        abstol, m, w, z, ldz, ifail);
}

template<typename T>
void dnevd(char jobz, char uplo, int n,
           T* a, int lda, Real<T>* w)
{
    LapackConnector::heevd(
        LapackConnector::ColMajor,
        jobz, uplo, n,
        a, lda, w);
}

template<typename T>
void potrf(char uplo, int n, T* A, int lda)
{
	LapackConnector::potrf(LapackConnector::ColMajor, uplo, n, A, lda);
}	


template<typename T>
void trtri(char uplo, char diag, int n, T* A, int lda)
{
    LapackConnector::trtri(LapackConnector::ColMajor, uplo, diag, n, A, lda);
}


template<typename T>
void getrf(int m, int n, T* A, int lda, int* ipiv)
{
    LapackConnector::getrf(LapackConnector::ColMajor, m, n, A, lda, ipiv);
}

template<typename T>
void getri(int n, T* A, int lda, int* ipiv)
{
    LapackConnector::getri(LapackConnector::ColMajor, n, A, lda, ipiv);
}

template<typename T>
void getrs(char trans, int n, int nrhs, T* A, int lda, int* ipiv, T* B, int ldb)
{
    LapackConnector::getrs(LapackConnector::ColMajor, trans, n, nrhs, A, lda, ipiv, B, ldb);
}

} // namespace lapackConnector
} // namespace container

#endif  // BASE_THIRD_PARTY_LAPACK_H_
