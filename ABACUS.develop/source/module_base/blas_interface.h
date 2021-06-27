#ifndef BLAS_INTERFACE_H
#define BLAS_INTERFACE_H

#include <complex>
using namespace std;

extern "C"
{
	// level 1: vector-vector operations	
	void dscal_(const int *n, const double *alpha, const double *x, const int *incx);
	void zscal_(int *n, complex<double> *alpha, complex<double>  *x, int *incx);

	void daxpy_(const int *n, const double* sa, const double* sx, const int* incx, const double* sy, const int* incy);
	void zaxpy_(const int *n, const complex<double> *sa, const complex<double> *sx, const int* incx, const complex<double> *sy, const int* incy);
	
	void dcopy_(int const *n, double *a, int const *incx, double *b, int const *incy); 

	void zdotc_(complex<double> *result, const int *n, const complex<double> *zx, 
	const int *incx, const complex<double> *zy, const int *incy);

	double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);

	// level 2: matrix-vector operations
	void dgemv_(const char *transa, const int *m, const int *n, const double *alpha,  const double *a,  
		const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
		
	void zgemv_(const char *trans, const int *m, const int *n, const complex<double> *alpha,
			const complex<double> *a, const int *lda, const complex<double> *x, const int *incx,
			const complex<double> *beta, complex<double> *y, const int *incy);

	// level 3: matrix-matrix operations
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
        const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
        const double *beta, double *c, const int *ldc);

    void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
        const complex<double> *alpha, const complex<double> *a, const int *lda, const complex<double> *b,
        const int *ldb, const complex<double> *beta, complex<double> *c, const int *ldc);

    void dzgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
        const complex<double> *alpha, const double *a, const int *lda, const complex<double> *b,
        const int *ldb, const complex<double> *beta, complex<double> *c, const int *ldc);

	void dsymm_(char *side, char *uplo, int *m, int *n, double *alpha,
	            double *a,  int *lda,  double *b, int *ldb, double *beta, double *c, int *ldc);

	void zhemm_(char *side, char *uplo, int *m, int *n,complex<double> *alpha,
		complex<double> *a,  int *lda,  complex<double> *b, int *ldb, complex<double> *beta, complex<double> *c, int *ldc);

	void dtrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
		double* alpha, double* a, int *lda, double*b, int *ldb);

	void ztrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
	complex<double>* alpha, complex<double>* a, int *lda, complex<double>*b, int *ldb);

	// computes all eigenvalues of a symmetric tridiagonal matrix
	// using the Pal-Walker-Kahan variant of the QL or QR algorithm.
	void dsterf_(int *n, double *d, double *e, int *info);

	void dstein_(int *n, double* d, double *e, int *m, double *w,
		int* block, int* isplit, double* z, int *lda, double *work,
		int* iwork, int* ifail, int *info);
 	void zstein_(int *n, double* d, double *e, int *m, double *w,
        int* block, int* isplit, complex<double>* z, int *lda, double *work,
        int* iwork, int* ifail, int *info);

	// computes the Cholesky factorization of a real symmetric
	// positive definite matrix A.
	void dpotf2_(char *uplo, int *n, double *a, int *lda, int *info);
	void zpotf2_(char *uplo,int *n,complex<double> *a, int *lda, int *info);

	void dsygs2_(int *itype, char *uplo, int *n, double *a, int *lda, double *b, int *ldb, int *info);
	void zhegs2_(int *itype, char *uplo, int *n, complex<double> *a, int *lda, complex<double> *b, int *ldb, int *info);

	void dlacpy_(char *uplo, int *m, int *n, double* a, int *lda, double *b, int *ldb);
	void zlacpy_(char *uplo, int *m, int *n, complex<double>* a, int *lda, complex<double> *b, int *ldb);

	void dlarfg_(int *n, double *alpha, double *x, int *incx, double *tau);
	void zlarfg_(int *n, complex<double> *alpha, complex<double> *x, int *incx, complex<double> *tau);

	void dger_(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
	void zgerc_(int *m, int *n, complex<double> *alpha,complex<double> *x, int *incx, complex<double> *y, int *incy,complex<double> *a, int *lda);

};

#endif
