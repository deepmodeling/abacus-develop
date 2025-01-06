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

    void pdtran_(const int* m, const int* n,
        const double* alpha, const double* a, const int* ia, const int* ja, const  int* desca,
        const double* beta, double* c, const int* ic, const int* jc, const int* descc);

	void pztranu_(const int *m,const int*n,
        const std::complex<double>* alpha, const std::complex<double>* a, const int* ia, const int* ja, const int* desca,
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

	void pssygvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, float* A, const int* ia, const int* ja, const int*desca, float* B, const int* ib, const int* jb, const int*descb,
		const float* vl, const float* vu, const int* il, const int* iu,
		const float* abstol, int* m, int* nz, float* w, const float*orfac, float* Z, const int* iz, const int* jz, const int*descz,
		float* work, int* lwork, int*iwork, int*liwork, int* ifail, int*iclustr, float*gap, int* info);
		
	void pchegvx_(const int* itype, const char* jobz, const char* range, const char* uplo,
		const int* n, std::complex<float>* A, const int* ia, const int* ja, const int*desca, std::complex<float>* B, const int* ib, const int* jb, const int*descb,
		const float* vl, const float* vu, const int* il, const int* iu,
		const float* abstol, int* m, int* nz, float* w, const float*orfac, std::complex<float>* Z, const int* iz, const int* jz, const int*descz,
		std::complex<float>* work, int* lwork, float* rwork, int* lrwork, int*iwork, int*liwork, int* ifail, int*iclustr, float*gap, int* info);


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

	// Scalapack wrappers to copy 2D blocks of data
	// more info: 
	// https://netlib.org/scalapack/explore-html/da/db5/pigemr_8c.html
	// https://netlib.org/scalapack/explore-html/dd/dcd/pdgemr_8c.html
	// https://netlib.org/scalapack/explore-html/d5/dd4/pzgemr_8c.html
	// https://netlib.org/scalapack/explore-html/d5/deb/psgemr_8c.html
	// https://netlib.org/scalapack/explore-html/d4/dad/pcgemr_8c.html
	void Cpigemr2d (int m, int n, int *ptrmyblock, int ia, int ja, int *ma, int *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
	void Cpdgemr2d (int m, int n, double *ptrmyblock, int ia, int ja, int *ma, double *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
	void Cpzgemr2d (int m, int n, std::complex<double> *ptrmyblock, int ia, int ja, int *ma, std::complex<double> *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
	void Cpsgemr2d (int m, int n, float *ptrmyblock, int ia, int ja, int *ma, float *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
	void Cpcgemr2d (int m, int n, std::complex<float> *ptrmyblock, int ia, int ja, int *ma, std::complex<float> *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
}

	template <typename T>
	struct block2d_data_type
	{
		constexpr static bool value = std::is_same<T, double>::value || std::is_same<T, std::complex<double>>::value || std::is_same<T, float>::value || std::is_same<T, std::complex<float>>::value || std::is_same<T, int>::value;
	};


	/**
	 * Copies a 2D block of data from matrix A to matrix B using the Scalapack library.
	 * This function supports different data types: double, std::complex<double>, float, std::complex<float>, and int.
	 *
	 * @tparam T The data type of the matrices A and B.
	 * @param M The number of rows of matrix A.
	 * @param N The number of columns of matrix A.
	 * @param A Pointer to the source matrix A.
	 * @param IA The starting row index of the block in matrix A.
	 * @param JA The starting column index of the block in matrix A.
	 * @param DESCA Descriptor array for matrix A.
	 * @param B Pointer to the destination matrix B.
	 * @param IB The starting row index of the block in matrix B.
	 * @param JB The starting column index of the block in matrix B.
	 * @param DESCB Descriptor array for matrix B.
	 * @param ICTXT The context identifier.
	 */
	template <typename T>
	typename std::enable_if<block2d_data_type<T>::value,void>::type Cpxgemr2d(int M, int N, T *A, int IA, int JA, int *DESCA, T *B, int IB, int JB, int *DESCB, int ICTXT)
	{
		if (std::is_same<T,double>::value) Cpdgemr2d(M, N, reinterpret_cast<double*>(A),IA, JA, DESCA,reinterpret_cast<double*>(B),IB,JB, DESCB,ICTXT);
		if (std::is_same<T,std::complex<double>>::value) Cpzgemr2d(M, N, reinterpret_cast<std::complex<double>*>(A),IA, JA, DESCA,reinterpret_cast<std::complex<double>*>(B),IB,JB, DESCB,ICTXT);
		if (std::is_same<T,float>::value) Cpsgemr2d(M, N, reinterpret_cast<float*>(A),IA, JA, DESCA,reinterpret_cast<float*>(B),IB,JB, DESCB,ICTXT);
		if (std::is_same<T,std::complex<float>>::value) Cpcgemr2d(M, N, reinterpret_cast<std::complex<float>*>(A),IA, JA, DESCA,reinterpret_cast<std::complex<float>*>(B),IB,JB, DESCB,ICTXT);
		if (std::is_same<T,int>::value) Cpigemr2d(M, N, reinterpret_cast<int*>(A),IA, JA, DESCA,reinterpret_cast<int*>(B),IB,JB, DESCB,ICTXT);
	};
	

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
            const double alpha,
            const double* A, const int IA, const int JA, const int* DESCA,
            const double* B, const int IB, const int JB, const int* DESCB,
            const double beta,
            double* C, const int IC, const int JC, const int* DESCC)
    {
        pdgemm_(&transa, &transb, &M, &N, &K, &alpha, A, &IA, &JA, DESCA,
            B, &IB, &JB, DESCB, &beta, C, &IC, &JC, DESCC);
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
