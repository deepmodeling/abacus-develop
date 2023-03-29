// =============================================================================
//                          C++ Header File
// Project:         LapackConnector
// File:            LapackConnector.hpp
// Author:          sltk
// Comment:         LapackConnector provide the connector to the fortran Lapack routine.
// Warning:         
// Start time:      2007-03-08
// Last modified:   2008-08-12 ywcui : add zhegvx
// 					2008-08-13 mohan : find bug,test.
// 					2008-09-03 mohan : Add zgesv
// =============================================================================

#ifndef LAPACKCONNECTOR_HPP
#define LAPACKCONNECTOR_HPP

#include <new>

//#include <exception>
#include <stdexcept>

#include <complex>
using namespace std;

#include "complexmatrix.h"


extern "C"
{

	void zpotrf_( const char* uplo,
				const int* n,
				const complex<double> *A,
				const int* lda,
				const int* infp);

	void zpotri_( const char* uplo,
				const int* n,
				const complex<double> *A,
				const int* lda,
				const int* infp);

	//=======================================================================
	//ZGESV computes the solution to a complex system of linear equations
 	//     A * X = B,
	//  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
	//
	//  The LU decomposition with partial pivoting and row interchanges is
	//  used to factor A as
	//     A = P * L * U,
	//  where P is a permutation matrix, L is unit lower triangular, and U is
	//  upper triangular.  The factored form of A is then used to solve the
	//  system of equations A * X = B.
	//
	//  A = LDA * N
	//  B = LDB * NRHS (exit X = B)
	//========================================================================
	void zgesv_(const int* n,
				const int* nrhs,
				const complex<double> *A,
				const int *lda,
				const int *ipiv,
				const complex<double> *B,
				const int* ldb,
				int* info
				);

	/*zhegv computes all the eigenvalues, and optionally, the eigenvectors
	  of a complex generalized Hermitian-definite eigenproblem, of the form
	  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
	  Here A and B are assumed to be Hermitian and B is also
	  positive definite.
	*/

	/*zhegvx computes selected eigenvalues, and optionally, eigenvectors
	  of a complex generalized Hermitian-definite eigenproblem, of the form
	  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
	  B are assumed to be Hermitian and B is also positive definite.
	  Eigenvalues and eigenvectors can be selected by specifying either a
	  range of values or a range of indices for the desired eigenvalues.
	*/

	void zhegv_(const int* itype,
                const char* jobz,
                const char* uplo,
                const int* n,
                complex<double>* a,
                const int* lda,
                complex<double>* b,
                const int* ldb,
                double* w,
                complex<double>* work,
                int* lwork,
                double* rwork,
                int* info	);

	void zheev_(const char* jobz,
				const char* uplo,
				const int* n,
				complex<double> *a,
				const int* lda,
				double* w,
				complex<double >* work,
				const int* lwork,
				double* rwork,
				int* info	);

	//===================================================================================
	// zhegvx	:
	// itype	:	1:	A*x = (lambda)*B*x
	// 				2:	A*B*x = (lambda)*x
	// 				3:	B*A*x = (lambda)*x
	// jobz		:	'N':  Compute eigenvalues only
	// 				'V':  Compute eigenvalues and eigenvectors.
	// range	:	'A': all eigenvalues will be found.
	//				'V': all eigenvalues in the half-open interval (VL,VU] will be found.
	//				'I': the IL-th through IU-th eigenvalues will be found.
	// uplo		:	'U':  Upper triangles of A and B are stored
	// 				'L':  Lower triangles of A and B are stored
	// N		: 	The order of the matrices A and B.  N >= 0.
	// A		:	(input/output) COMPLEX*16 array, dimension (LDA, N)
    //            	On entry, the Hermitian matrix A.  If UPLO = 'U',  the  leading
    //           	N-by-N upper triangular part of A contains the upper triangular
    //           	part of the matrix A.  If UPLO = 'L', the leading N-by-N  lower
    //           	triangular  part of A contains the lower triangular part of the
    //          	matrix A.
	// LDA		:	The leading dimension of the array A.  LDA >= max(1,N).
	// LDB		:	On  entry,  the Hermitian matrix B.  If UPLO = 'U', the leading
    //            	N-by-N upper triangular part of B contains the upper triangular
    //           	part  of the matrix B.  If UPLO = 'L', the leading N-by-N lower
    //           	triangular part of B contains the lower triangular part of  the
    //           	matrix B.
	// vl,vu	:	the lower and upper bounds of the interval to be searched for 
	//				eigenvalues. VL < VU.  Not referenced if RANGE = 'A' or 'I'.
	// il,iu	:	the indices (in ascending order) of the smallest and largest 
	//				eigenvalues to be  returned. 1  <= IL <= IU <= N, if N > 0; 
	//				IL = 1 and IU = 0 if N = 0.  Not referenced if RANGE = 'A' or 'V'.
	// abstol	:	The absolute error tolerance for the eigenvalues
	// m		:	The  total number of eigenvalues found ,0 <= M <= N.  If RANGE
	//              = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
	// w		:	The first  M  elements  contain  the  selected  eigenvalues  in
	//              ascending order
	// z		:   If  JOBZ  = 'N', then Z is not referenced.  If JOBZ = 'V', then
    //           	if INFO = 0, the first M columns of Z contain  the  orthonormal
    //          	eigenvectors  of  the  matrix  A  corresponding to the selected
    //          	eigenvalues, with the i-th column of Z holding the  eigenvector
    //           	associated  with W(i).  The eigenvectors are normalized as fol-
    //           	lows: if ITYPE  =  1  or  2,  Z**T*B*Z  =  I;  if  ITYPE  =  3,
    //           	Z**T*inv(B)*Z = I.
	//				If  an  eigenvector  fails  to  converge, then that column of Z
    //		        contains the latest approximation to the eigenvector,  and  the
    //	            index  of the eigenvector is returned in IFAIL.  Note: the user
    //              must ensure that at least max(1,M) columns are supplied in  the
    //              array  Z;  if RANGE = 'V', the exact value of M is not known in
    //              advance and an upper bound must be used 
	// LDZ		:	The leading dimension of the array Z.  LDZ >= 1, and if JOBZ  =
    //              'V', LDZ >= max(1,N).
	// work		:	On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	// lwork	:	LWORK >= max(1,2*N). For optimal efficiency, LWORK >= (NB+1)*N, 
	//				where NB is the blocksize for ZHETRD returned by ILAENV.
	// rwork	:   dimension (7*N)
	// iwork	:	dimension (5*N)
	// info		:	= 0:  successful exit
	// 				< 0:  if INFO = -i, the i-th argument had an illegal value
	// 				> 0:  ZPOTRF or ZHEEVX returned an error code:
	//===================================================================================

	/*
	void zhegvx_(const int* itype,
				 const char* jobz,
				 const char* range,	
				 const char* uplo,
				 const int* N,
				 complex<double> *A,
				 const int* LDA,
				 complex<double> *B,
				 const int* LDB,
				 const double* vl,
				 const double* vu,
				 const int* il,
				 const int* iu,
				 const double* abstol,
				 const int* m,
				 double* w,
				 complex<double> *z,
				 const int *LDZ,
				 complex<double> *work,
				 const int* lwork,
				 double* rwork,
				 int* iwork,
				 int* ifail,
				 int* info
				 );
	*/

	
	// Peize Lin add 2015-11-20	
	void daxpy_(int *N, double *alpha, double *X, int *incX, double *Y, int *incY);
	void dscal_(int *N, double *alpha, double *X, int *incX);

	// Peize Lin add 2015-12-07
	void zgemm_(const char *TransA, const char *TransB,
                       const int *M, const int *N, const int *K,
					   const double *alpha, const void *A, const int *lda, const void *B, const int *ldb,
                       const double *beta, void *C, const int *ldc);
	void zsymm_(const char *Side, const char *Uplo,
                       const int *M, const int *N,
					   const double *alpha, const void *A, const int *lda, const void *B, const int *ldb,
                       const double *beta, void *C, const int *ldc);					   
	void zhemm_(const char *Side, const char *Uplo,
                       const int *M, const int *N,
					   const double *alpha, const void *A, const int *lda, const void *B, const int *ldb,
                       const double *beta, void *C, const int *ldc);					   
	
};

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
class LapackConnector
{
private:
	// Transpose the complex matrix to the fortran-form real-complex array.
	static inline
	complex<double>* transpose(const ComplexMatrix& a, const int n, const int lda)
	{	
		complex<double>* aux = new complex<double>[lda*n];
		for(int i = 0; i < n; ++i)
		{	for(int j = 0; j < lda; ++j)
			{	
				aux[i*lda+j] = a(j,i);		// aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
			}
		}
		return aux;
	}

	// Transpose the fortran-form real-complex array to the complex matrix.
	static inline
	void transpose(const complex<double>* aux, ComplexMatrix& a, const int n, const int lda)
	{	
		for(int i = 0; i < n; ++i)
		{	for(int j = 0; j < lda; ++j)
			{	
				a(j, i) = aux[i*lda+j];		// aux[i*lda+j] means aux[i][j] in semantic, not in syntax!
			}
		}
	}

	static inline
	char change_side(const char &side)
	{
		switch(side)
		{
			case 'L': return 'R';
			case 'R': return 'L';
			default: throw invalid_argument("Side must be 'L' or 'R'");
		}		
	}
	
	static inline
	char change_uplo(const char &uplo)
	{
		switch(uplo)
		{
			case 'U': return 'L';
			case 'L': return 'U';
			default: throw invalid_argument("Uplo must be 'U' or 'L'");
		}		
	}


	
public:
	static inline
	void zgesv (const int n,
				const int nrhs,
				ComplexMatrix &A,
				const int lda,
				const int *ipiv,
				complex<double> *B,
				const int ldb,
				int* info
				)
	{
		complex<double> *aux = LapackConnector::transpose(A, n, lda);

		zgesv_(&n, &nrhs, aux, &lda, ipiv, B, &ldb, info);

		LapackConnector::transpose(aux, A, n, lda);
	}

	static inline
	void zpotrf(char uplo,
				int n,
				ComplexMatrix &a,
				const int lda,
				int& info)
	{
/*		complex<double> *aux = LapackConnector::transpose(a, n, lda);
		zpotrf_( &uplo, &n, aux, &lda, &info);
		LapackConnector::transpose(aux, a, n, lda);
		delete[] aux;
*/		// Peize Lin update 2015-12-27
		const char uplo_changed(change_uplo(uplo));
		zpotrf_( &uplo_changed, &n, a.c, &lda, &info);
	}

	static inline
	void zpotri(char uplo,
				int n,
				ComplexMatrix &a,
				const int lda,
				int& info)
	{
/*		complex<double> *aux = LapackConnector::transpose(a, n, lda);
		zpotri_( &uplo, &n, aux, &lda, &info);
		LapackConnector::transpose(aux, a, n, lda);
		delete[] aux;
*/		// Peize Lin update 2015-12-27
		const char uplo_changed(change_uplo(uplo));
		zpotri_( &uplo_changed, &n, a.c, &lda, &info);		
	}


	// wrap function of fortran lapack routine zhegv.
    static inline
	void zhegv(	const int itype,
                const char jobz,
                const char uplo,
                const int n,
                ComplexMatrix& a,
                const int lda,
                ComplexMatrix& b,
                const int ldb,
                double* w,
                complex<double>* work,
                int lwork,
                double* rwork,
                int& info	)
	{	// Transpose the complex matrix to the fortran-form real-complex array.
		complex<double>* aux = LapackConnector::transpose(a, n, lda);
		complex<double>* bux = LapackConnector::transpose(b, n, ldb);
		
		// call the fortran routine
		zhegv_(&itype, &jobz, &uplo, &n, aux, &lda, bux, &ldb, w, work, &lwork, rwork, &info);

		// Transpose the fortran-form real-complex array to the complex matrix.
		LapackConnector::transpose(aux, a, n, lda);
		LapackConnector::transpose(bux, b, n, ldb);

		// free the memory.
		delete[] aux;
		delete[] bux;
	}

	// wrap function of fortran lapack routine zheev.
	static inline
	void zheev( const char jobz,
				const char uplo,
				const int n,
				ComplexMatrix& a,
				const int lda,
				double* w,
				complex< double >* work,
				const int lwork,
				double* rwork,
				int* info	)
	{	// Transpose the complex matrix to the fortran-form real-complex array.
		complex<double> *aux = LapackConnector::transpose(a, n, lda);

		// call the fortran routine
		zheev_(&jobz, &uplo, &n, aux, &lda, w, work, &lwork, rwork, info);

		// Transpose the fortran-form real-complex array to the complex matrix.
		LapackConnector::transpose(aux, a, n, lda);

		// free the memory.
		delete[] aux;
	}
	
	// wrap function of fortran lapack routine zheev.	
	static inline
	void zhegvx( const int itype,
				 const char jobz,
				 const char range,
				 const char uplo,
				 const int n,
				 ComplexMatrix& a,
				 const int lda,
				 ComplexMatrix& b,
				 const int ldb,
				 const double vl,
				 const double vu,
				 const int il,
				 const int iu,
				 const double abstol,
				 const int m,
				 double* w,
				 ComplexMatrix& z,
				 const int ldz,
				 complex<double>* work,
				 const int lwork,
				 double* rwork,
				 int* iwork,
				 int* ifail,
				 int& info)
	{
		 // Transpose the complex matrix to the fortran-form real-complex array.
		complex<double>* aux = LapackConnector::transpose(a, n, lda);
		complex<double>* bux = LapackConnector::transpose(b, n, ldb);
		complex<double>* zux = LapackConnector::transpose(z, n, ldz);
		
		// call the fortran routine
//	    zhegvx_(&itype, &jobz, &range, &uplo, &n, aux, &lda, bux, &ldb, &vl,
//	    &vu, &il,&iu, &abstol, &m, w, zux, &ldz, work, &lwork, rwork, iwork, ifail, &info);
	    
	    // Transpose the fortran-form real-complex array to the complex matrix
		LapackConnector::transpose(zux, z, n, ldz);	 
           
		// free the memory.
		delete[] aux;
		delete[] bux;   
 		delete[] zux;
		
     }             

	
	// Peize Lin add 2015-11-20
	static inline 
	void daxpy( int N,  double alpha,  double *X,  int incX, double *Y,  int incY)
	{
		daxpy_(&N, &alpha, X, &incX, Y, &incY);
	}	
	
	static inline 
	void dscal( int N,  double alpha, double *X,  int incX)
	{
		dscal_(&N, &alpha, X, &incX);
	}

	// Peize Lin add 2015-12-07
	// multiply two matrix. C=A*B
	static inline 
	void zgemm(const char TransA, const char TransB,
			   const int M, const int N, const int K,
			   const double alpha, const void *A, const int lda, const void *B, const int ldb,
			   const double beta, void *C, const int ldc)				   
	{
		zgemm_(&TransA, &TransB,
			   &N, &M, &K,
			   &alpha, B, &ldb, A, &lda,
			   &beta, C, &ldc);
	}

	// multiply two matrix. C=A*B('L'eft) or C=B*A('R'ight). A^T=A, stored in 'U'pper or 'L'ower
	static inline 
	void zsymm(const char side, const char uplo,
			   const int M, const int N,
			   const double alpha, const void *A, const int lda, const void *B, const int ldb,
			   const double beta, void *C, const int ldc)				   
	{
		const char &side_changed(change_side(side));
		const char &uplo_changed(change_uplo(uplo));
		zsymm_(&side_changed, &uplo_changed,
			   &N, &M,
			   &alpha, A, &lda, B, &ldb,
			   &beta, C, &ldc);
	}
	
	// multiply two matrix. C=A*B('L'eft) or C=B*A('R'ight). A^H=A, stored in 'U'pper or 'L'ower
	static inline 
	void zhemm(const char side, const char uplo,
			   const int M, const int N,
			   const double alpha, const void *A, const int lda, const void *B, const int ldb,
			   const double beta, void *C, const int ldc)				   
	{
		const char &side_changed(change_side(side));
		const char &uplo_changed(change_uplo(uplo));
		zhemm_(&side_changed, &uplo_changed,
			   &N, &M,
			   &alpha, A, &lda, B, &ldb,
			   &beta, C, &ldc);
	}
};

#endif  // LAPACKCONNECTOR_HPP
