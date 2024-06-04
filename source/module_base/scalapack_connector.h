#ifndef SCALAPACK_CONNECTOR_H
#define SCALAPACK_CONNECTOR_H

#ifdef __MPI

#include <complex>

extern "C"
{
	int numroc_( const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs );
	void descinit_( 
		int *desc, 
		const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, 
		const int *ictxt, const int *lld, int *info);

	void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info);
//	void pzpotrf_(char *uplo, int *n, double _Complex *a, int *ia, int *ja, int *desca, int *info);
	void pzpotrf_(char *uplo, int *n, std::complex<double> *a, int *ia, int *ja, int *desca, int *info);

	void pdtran_(int *m , int *n ,
		double *alpha , double *a , int *ia , int *ja , int *desca ,
		double *beta ,  double *c , int *ic , int *jc , int *descc );

	void pztranu_(const int *m,const int*n,
		const std::complex<double> *alpha , std::complex<double> *a , const int *ia , const int *ja ,const int *desca ,
		const std::complex<double> *beta ,  std::complex<double> *c , const int *ic ,const int *jc ,const int *descc);

	void pzgemv_(
		const char *transa,
		const int *M, const int *N,
		const double *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB, const int *K, 
		const double *beta, std::complex<double> *C, const int *IC, const int *JC, const int *DESCC,const int *L);
	void pdgemv_(
		const char *transa,
		const int *M, const int *N,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB, const int *K, 
		const double *beta, double *C, const int *IC, const int *JC, const int *DESCC,const int *L);
	// C = a * A.? * B.? + b * C
	void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);
	void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
	void pdsymm_(char *side , char *uplo , int *m , int *n ,
		double *alpha , double *a , int *ia , int *ja , int *desca ,
		double *b , int *ib , int *jb , int *descb ,
		double *beta ,  double *c , int *ic , int *jc , int *descc );
	void pdtrmm_(char *side , char *uplo , char *transa , char *diag , int *m , int *n ,
		double *alpha , double *a , int *ia , int *ja , int *desca ,
		double *b , int *ib , int *jb , int *descb );
//	void pztrmm_(char *side , char *uplo , char *transa , char *diag , int *m , int *n ,
//		double *alpha , double _Complex *a , int *ia , int *ja , int *desca ,
//		double _Complex *b , int *ib , int *jb , int *descb );
	void pztrmm_(char *side , char *uplo , char *transa , char *diag , int *m , int *n ,
		std::complex<double> *alpha , std::complex<double> *a , int *ia , int *ja , int *desca ,
		std::complex<double> *b , int *ib , int *jb , int *descb );

	void pzgetrf_(
		const int *M, const int *N, 
		std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		int *ipiv,  int *info);

	void pdsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, double* A, const int* ia, const int* ja, const int*desca, double* B, const int* ib, const int* jb, const int*descb,
		const double* vl, const double* vu, const int* il, const int* iu,
		const double* abstol, int* m, int* nz, double* w, const double*orfac, double* Z, const int* iz, const int* jz, const int*descz,
		double* work, int* lwork, int*iwork, int*liwork, int* ifail, int*iclustr, double*gap, int* info);
	void pzhegvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, std::complex<double>* A, const int* ia, const int* ja, const int*desca, std::complex<double>* B, const int* ib, const int* jb, const int*descb,
		const double* vl, const double* vu, const int* il, const int* iu,
		const double* abstol, int* m, int* nz, double* w, const double*orfac, std::complex<double>* Z, const int* iz, const int* jz, const int*descz,
		std::complex<double>* work, int* lwork, double* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, double*gap, int* info);

	void pzgetri_(
		const int *n, 
		const std::complex<double> *A, const int *ia, const int *ja, const int *desca,
		int *ipiv, const std::complex<double> *work, const int *lwork, const int *iwork, const int *liwork, const int *info);

    void pzgeadd_(
		const char *transa,
		const int *m, const int *n,
		const std::complex<double> *alpha,
		const std::complex<double> *a, const int *ia, const int *ja, const int *desca,
                		const std::complex<double> *beta,
		const std::complex<double> *c, const int *ic, const int *jc, const int *descc);

    void pztranc_(
		const int *M, const int *N,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
	
    void pdgemr2d_(const int *M, const int *N,
	    double *A, const int *IA, const int *JA, const int *DESCA, 
		double *B, const int *IB, const int *JB, const int *DESCB,
		const int *ICTXT);			   
		
    void pzgemr2d_(const int *M, const int *N,
	    std::complex<double> *A, const int *IA, const int *JA, const int *DESCA, 
		std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const int *ICTXT);
}

class ScalapackConnector
{
public:
	static inline
	void geadd(
		const char transa,
		const int m, const int n,
		const std::complex<double> alpha,
		const std::complex<double> *a, const int ia, const int ja, const int *desca,
		const std::complex<double> beta,
		const std::complex<double> *c, const int ic, const int jc, const int *descc)
	{
		pzgeadd_(&transa, &m, &n, &alpha, a, &ia, &ja, desca, &beta, c, &ic, &jc, descc);
	}

	static inline
	void gemm(
		const char transa, const char transb,
		const int M, const int N, const int K,
		const std::complex<double> alpha,
		const std::complex<double> *A, const int IA, const int JA, const int *DESCA,
		const std::complex<double> *B, const int IB, const int JB, const int *DESCB,
		const std::complex<double> beta,
		std::complex<double> *C, const int IC, const int JC, const int *DESCC)
	{
		pzgemm_(&transa, &transb, &M, &N, &K, &alpha, A, &IA, &JA, DESCA, 
			B, &IB, &JB, DESCB, &beta, C, &IC, &JC, DESCC);
	}

	static inline
	void getrf(
		const int M, const int N, 
		std::complex<double> *A, const int IA, const int JA, const int *DESCA,
		int *ipiv,  int *info) //fix a bug: info is output and we must use int*
	{
		pzgetrf_(&M, &N, A, &IA, &JA, DESCA, ipiv, info);
	}

	static inline
	void getri(
		const int n, 
		const std::complex<double> *A, const int ia, const int ja, const int *desca, int *ipiv, 
		const std::complex<double> *work, const int *lwork, const int *iwork, const int *liwork, int *info)
	{
		pzgetri_(&n, A, &ia, &ja, desca, ipiv, work, lwork, iwork, liwork, info);
	}

	static inline
	void tranu(
		const int m, const int n,
		const std::complex<double> alpha , std::complex<double> *a , const int ia , const int ja , const int *desca,
		const std::complex<double> beta ,  std::complex<double> *c , const int ic , const int jc , const int *descc)
	{
		pztranu_(&m, &n, &alpha, a, &ia, &ja, desca, &beta, c, &ic, &jc, descc);
	}
};

#endif // __MPI

#endif
