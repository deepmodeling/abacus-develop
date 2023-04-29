#ifndef CONTAINER_KERNELS_THIRD_PARTY_LAPACK_CONNECTOR_H_
#define CONTAINER_KERNELS_THIRD_PARTY_LAPACK_CONNECTOR_H_

#include <complex>

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
int ilaenv_(int* ispec,const char* name,const char* opts,
            const int* n1,const int* n2,const int* n3,const int* n4);


// solve the generalized eigenproblem Ax=eBx, where A is Hermitian and complex couble
// zhegv_ & zhegvd_ returns all eigenvalues while zhegvx_ returns selected ones

void chegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda,
             const std::complex<float>* b, const int* ldb, float* w,
             std::complex<float>* work, int* lwork, float* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

void zhegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<double>* a, const int* lda,
             const std::complex<double>* b, const int* ldb, double* w,
             std::complex<double>* work, int* lwork, double* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

void cheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<float> *a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu, const float* abstol,
             const int* m, float* w, std::complex<float> *z, const int *ldz,
             std::complex<float> *work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);

void zheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<double> *a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu, const double* abstol,
             const int* m, double* w, std::complex<double> *z, const int *ldz,
             std::complex<double> *work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);
}

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
class LapackConnector
{
private:

    // wrap function of fortran lapack routine zhegvd.
    static inline
    void zhegvd(const int itype, const char jobz, const char uplo, const int n,
                std::complex<double>* a, const int lda,
                const std::complex<double>* b, const int ldb, double* w,
                std::complex<double>* work, int lwork, double* rwork, int lrwork,
                int* iwork, int liwork, int info)
    {
        zhegvd_(&itype, &jobz, &uplo, &n,
                a, &lda, b, &ldb, w,
                work, &lwork, rwork, &lrwork,
                iwork, &liwork, &info);
    }


    // wrap function of fortran lapack routine zheevx.
    static inline
    void zheevx( const int itype, const char jobz, const char range, const char uplo, const int n,
                 std::complex<double>* a, const int lda,
                 const double vl, const double vu, const int il, const int iu, const double abstol,
                 const int m, double* w, std::complex<double>* z, const int ldz,
                 std::complex<double>* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
    {
        zheevx_(&jobz, &range, &uplo, &n,
                a, &lda, &vl, &vu, &il, &iu,
                &abstol, &m, w, z, &ldz,
                work, &lwork, rwork, iwork, ifail, &info);
    }

public:
// wrap function of fortran lapack routine zhegvd. (pointer version)
static inline
void xhegvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<float>* a, const int lda,
            const std::complex<float>* b, const int ldb, float* w,
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
void xhegvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<double>* a, const int lda,
            const std::complex<double>* b, const int ldb, double* w,
            std::complex<double>* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    zhegvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
}

    static inline
    int ilaenv( int ispec, const char *name,const char *opts,const int n1,const int n2,
                const int n3,const int n4)
    {
        const int nb = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
        return nb;
    }

// wrap function of fortran lapack routine zheevx.
static inline
void xheevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             std::complex<float>* a, const int lda,
             const float vl, const float vu, const int il, const int iu, const float abstol,
             const int m, float* w, std::complex<float>* z, const int ldz,
             std::complex<float>* work, const int lwork, float* rwork, int* iwork, int* ifail, int info)
{
    cheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}

// wrap function of fortran lapack routine zheevx.
static inline
void xheevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             std::complex<double>* a, const int lda,
             const double vl, const double vu, const int il, const int iu, const double abstol,
             const int m, double* w, std::complex<double>* z, const int ldz,
             std::complex<double>* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
{
    zheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}
};
#endif  // CONTAINER_KERNELS_THIRD_PARTY_LAPACK_CONNECTOR_H_
