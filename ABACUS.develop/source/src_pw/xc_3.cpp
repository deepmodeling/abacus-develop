/* myfunc.cpp */
// from LPACK
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

using namespace std;
#include "myfunc.h"
#include "module_base/blas_connector.h"
#include "../src_pw/global.h"

// dscal compute y = alpha * y, where alpha is a scalar and
void dscal(const int n,const double &alpha,double *y,const int incy)
{
    if (incy <= 0 || n < 1)
    {
        cout << "\n error in dscal,incy <= 0 or n < 1";
        return;
    }

    for (int i = 0; i < n; i += incy)
    {
        y[i] *= alpha;
    }

    return;
}

// a(i,:) = alpha * a(i,:) where a is a matrix
// i line
void dscal(const double &alpha,matrix &a,const int i)
{
    int nc = a.nc;
    int nr = a.nr;

    if (nc <= 0 || nr <= 0 || i < 0 || i >= nr)
    {
        cerr << "\n error in dscal,nc <= 0 or nr <= 0 or i < 0 or i >= nr, ";
        return;
    }

    for (int j = 0;j < nc;j++)
    {
        a(i,j) = a(i,j) * alpha;
    }
}

// daxpy compute y := alpha * x + y where alpha is a scalar and x and y
void daxpy(const int n, const double &alpha, const double *x, const int incx, double *y, const int incy)
{
    if (n < 1 || incy <= 0 || incx <= 0)
    {
        cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }
    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] += alpha * x[ix];
    }
    return;
}

void zaxpy(int n, double alpha, complex < double> *x, int incx, complex < double> *y, int incy)
{

    if (n < 1 || incy <= 0 || incx <= 0)
    {
        cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)  // y := alpha * x + y
    {
        y[iy] += alpha * x[ix];
    }
}


//-------------------------------------------------------------
void zaxpy(int n, complex < double> alpha,  complex < double> *x,
           int incx, complex < double> *y, int incy)
{
    // zaxpy compute y := alpha * x + y where alpha is a scalar and
    // x and y are n-vectors.
    // DOUBLE COMPLEX ALPHA
    // DOUBLE COMPLEX X(*), Y(*)
    // INTEGER N, INCX, INCY

    // ARGUMENTS

    // N (input)
    //		On entry, N specifies the number  of  elements  in
    //		the  vector.   N must be at least one for the sub-
    //		routine to have any visible effect.  Unchanged  on
    //		exit.

    // ALPHA (input)
    //		On  entry,  ALPHA  specifies  the  scalar   alpha.
    //		Unchanged on exit.

    // X (input)
    //		array of DIMENSION at least ( 1 + ( n -  1  )*abs(
    //		INCX  )  ).  Before entry, the incremented array X
    //		must contain the vector x.  Unchanged on exit.

    // INCX (input)
    //		On entry, INCX specifies  the  increment  for  the
    //		elements of X. Unchanged on exit.

    // Y (input/output)
    //		array of DIMENSION at least ( 1 + ( n -  1  )*abs(
    //		INCY  ) ).  On entry, the incremented array Y must
    //		contain the vector y. On exit, Y is overwritten by
    //		the updated vector y.

    // INCY (input)
    //		On entry, INCY specifies  the  increment  for  the
    //		elements of Y. Unchanged on exit.

    if (n < 1 || incy <= 0 || incx <= 0)
    {
        cerr << "\n error in daxpy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)  // y := alpha * x + y
    {
        y[iy] += alpha * x[ix];
    }

} // end zaxpy

// y(i,:) = alpha * x + y(i,:) where y is a matrix
void zaxpy(double alpha,
           complex < double> *x,
           ComplexMatrix &y,
           int i)
{
    int nr, nc;
    nr = y.nr;
    nc = y.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y(i, j) += alpha * x[j];
    }

//    cout << "\n End daxpy() " << endl;
}

// y = alpha * x(i,:) + y
void zaxpy(double alpha,
           const ComplexMatrix &x,
           int i,
           complex < double> *y)
{
    int nr, nc;
    nr = x.nr;
    nc = x.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y[j] += alpha * x(i, j);
    }
}

void zaxpy(complex < double> alpha,
           const ComplexMatrix &x,
           int i,
           complex < double> *y)
{
    int nr, nc;
    nr = x.nr;
    nc = x.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in daxpy, nr < 1 or nc < 1 or i out of range, ";
        return;
    }

    for (int j = 0;j < nc;j ++) // y := alpha * x + y
    {
        y[j] += alpha * x(i, j);
    }
}

// y(j,:) = alpha * x(i,:) + y(j,:)
void zaxpy(complex < double> alpha,
           const ComplexMatrix &x,
           int i,
           ComplexMatrix &y,
           int j)
{
    int nr, nc;
    nr = y.nr;
    nc = y.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr || j < 0 || j >= nr)
    {
        cerr << "\n error in daxpy, nr < 1 or nc < 1 or i or j out of range, ";
        return;
    }

    for (int k = 0;k < nc;k ++) // y := alpha * x + y
    {
        y(j, k) += alpha * x(i, k);
    }
}

//-----------------------------------------------------------------
void dcopy(int n, double *x, int incx, double *y, int incy)
{
    // dcopy Copy x to y where x and y are n-vectors.
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }
} // end dcopy

//------------------------------------------------------------------
void dcopy(int n, complex < double> *x, int incx, complex < double> *y, int incy)
{
    // zcopy Copy x to y where x and y are n-vectors.
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }

} // end dcopy

void dcopy(int n, int *x, int incx, int *y, int incy)
{
    if (n < 1 || incx <= 0 || incy <= 0)
    {
        cerr << "\n error in dcopy, n < 1 or incx <= 0 or incy <= 0, ";
        return;
    }

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        y[iy] = x[ix];
    }
}

/* Copy a(i,:) to y where x is matrix, and y are n-vectors. */
void dcopy(const matrix &a,
           int i,
           double *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy((matrix a, int i, double *y)), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = a(i, j);
    }

} // end dcopy

//-------------------------------------
void dcopy(const matrix &a,
           int i,
           int *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(matrix a, int i, int *y),"
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = (int) a(i, j);
    }

} // end dcopy

//-----------------------------------------------------------------
void dcopy(const ComplexMatrix &a,
           int i,
           complex < double> *y)
{
    // dcopy Copy a(i,:) to y where a is complex matrix, and y are n-vectors.
    const int nr = a.nr;
    const int nc = a.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(ComplexMatrix a, int i, complex < double> *),"
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        y[j] = a(i, j);
    }
} // end dcopy

// ------------------------------------
void dcopy(double *x, matrix &b, int i)
{
    // copy x to ith row of b where b is a matrix and x is a vector
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(double *x, matrix &b, int i), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        b(i, j) = x[j];
    }
}

void dcopy(complex < double> *x, ComplexMatrix &b, int i)
{
    // copy x to ith row of b where b is a complex matrix and x is a vector
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(double> *x, ComplexMatrix &b, int i), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int j = 0;j < nc; j++)
    {
        b(i, j) = x[j];
    }
}

// b(j,:) = a(i,:)
void dcopy(const matrix &a,
           int i,
           matrix &b,
           int j)
{
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(matrix a, int i,matrix &b, int ), "
             << "nr or nc < 1 or i out of range ";
        return;
    }

    for (int k = 0;k < nc; k++)
    {
        b(j, k) = a(i, k);
    }
}

void dcopy(const ComplexMatrix &a,
           int i,
           ComplexMatrix &b,
           int j)
{
    int nr, nc;
    nr = b.nr;
    nc = b.nc;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dcopy(ComplexMatrix a, int i,ComplexMatrix &b, int j), "
             << " nr or nc < 1 or i out of range ";
        return;
    }

    for (int k = 0;k < nc; k++)
    {
        b(j, k) = a(i, k);
    }
}

void dcopy(int n,
           const ComplexMatrix &a,
           int inca,
           ComplexMatrix &b,
           int i,
           int incb)
{
    cout << "\n do nothing, in dcopy() ";
}



//-------------------------------------------------------------------
//void dsytrf(char uplo,  int  n,  double  *a,  int  lda,  int *ipivot,
//			int *info)
void dsytrf(char , int iter_used, matrix beta, int maxter, int *iwork,
            double *work, int , int &info)
{
    // dsytrf computes the factorization of a real symmetric matrix
    // A  using  the  Bunch-Kaufman  diagonal pivoting method.  The
    // form of the factorization is
    //		A = U*D*U**T  or  A = L*D*L**T
    // where U (or L) is a product of permutation  and  unit  upper
    // (lower)  triangular  matrices,  and D is symmetric and block
    // diagonal with 1-by-1 and 2-by-2 diagonal blocks.
    // This is the blocked version of the algorithm, calling  Level
    // 3 BLAS.

    // CHARACTER * 1 UPLO
    // INTEGER N, LDA, LDWORK, INFO
    // INTEGER IPIVOT(*)
    // DOUBLE PRECISION A(LDA,*), WORK(*)
    /***********************************************************
    // UPLO (input)
    //		= 'U':  Upper triangle of A is stored;
    //		= 'L':  Lower triangle of A is stored.

    // N (input) The order of the matrix A.  N >= 0.

    // A (input/output)
               On entry, the symmetric matrix A.  If UPLO =  'U',
               the leading N-by-N upper triangular part of A con-
               tains the upper triangular part of the  matrix  A,
               and the strictly lower triangular part of A is not
               referenced.  If UPLO =  'L',  the  leading  N-by-N
               lower triangular part of A contains the lower tri-
               angular part of the matrix  A,  and  the  strictly
               upper triangular part of A is not referenced.

               On exit, the block diagonal matrix D and the  mul-
               tipliers  used  to  obtain  the factor U or L (see
               below for further details).

     LDA (input)
               The leading dimension of  the  array  A.   LDA  >=
               max(1,N).

     IPIVOT (output)
               Details of the interchanges and the  block  struc-
               ture  of  D.   If  IPIVOT(k)  >  0,  then rows and
               columns k  and  IPIVOT(k)  were  interchanged  and
               D(k,k)  is a 1-by-1 diagonal block.  If UPLO = 'U'
               and IPIVOT(k) = IPIVOT(k-1) <  0,  then  rows  and
               columns  k-1  and -IPIVOT(k) were interchanged and
               D(k-1:k,k-1:k) is a  2-by-2  diagonal  block.   If
               UPLO  =  'L' and IPIVOT(k) = IPIVOT(k+1) < 0, then
               rows and columns k+1 and  -IPIVOT(k)  were  inter-
               changed  and  D(k:k+1,k:k+1)  is a 2-by-2 diagonal
               block.

     WORK (workspace)
               On exit, if INFO = 0, WORK(1) returns the  optimal
               LDWORK.

     LDWORK (input)
               The length of WORK.  LDWORK >=1.  For best perfor-
               mance  LDWORK  >= N*NB, where NB is the block size
               returned by ILAENV.

               If LDWORK = -1, then a workspace query is assumed;
               the  routine  only  calculates the optimal size of
               the WORK array, returns this value  as  the  first
               entry  of  the  WORK  array,  and no error message
               related to LDWORK is issued by XERBLA.

     INFO (output)
               = 0:  successful exit
               < 0:  if INFO = -i, the i-th argument had an ille-
               gal value
               > 0:  if INFO = i, D(i,i) is  exactly  zero.   The
               factorization  has  been  completed, but the block
               diagonal matrix D is exactly singular,  and  divi-
               sion  by  zero will occur if it is used to solve a
               system of equations.

    FURTHER DETAILS

     If UPLO = 'U', then A = U*D*U', where
        U = P(n)*U(n)* ... *P(k)U(k)* ...,
     i.e., U is a product of terms P(k)*U(k), where  k  decreases
     from  n  to  1 in steps of 1 or 2, and D is a block diagonal
     matrix with 1-by-1 and 2-by-2 diagonal blocks D(k).  P(k) is
     a  permutation matrix as defined by IPIVOT(k), and U(k) is a
     unit upper triangular matrix,  such  that  if  the  diagonal
     block D(k) is of order s (s = 1 or 2), then

                (   I    v    0   )   k-s
        U(k) =  (   0    I    0   )   s
                (   0    0    I   )   n-k
                   k-s   s   n-k

     If s = 1, D(k) overwrites A(k,k), and  v  overwrites  A(1:k-
     1,k).   If s = 2, the upper triangle of D(k) overwrites A(k-
     1,k-1), A(k-1,k), and A(k,k), and  v  overwrites  A(1:k-2,k-
     1:k).

     If UPLO = 'L', then A = L*D*L', where
        L = P(1)*L(1)* ... *P(k)*L(k)* ...,
     i.e., L is a product of terms P(k)*L(k), where  k  increases
     from  1  to  n in steps of 1 or 2, and D is a block diagonal
     matrix with 1-by-1 and 2-by-2 diagonal blocks D(k).  P(k) is
     a  permutation matrix as defined by IPIVOT(k), and L(k) is a
     unit lower triangular matrix,  such  that  if  the  diagonal
     block D(k) is of order s (s = 1 or 2), then

                (   I    0     0   )  k-1
        L(k) =  (   0    I     0   )  s
                (   0    v     I   )  n-k-s+1
                   k-1   s  n-k-s+1

     If  s  =  1,  D(k)  overwrites  A(k,k),  and  v   overwrites
     A(k+1:n,k).  If s = 2, the lower triangle of D(k) overwrites
     A(k,k),  A(k+1,k),  and   A(k+1,k+1),   and   v   overwrites
     A(k+2:n,k:k+1).
    ******************************************************************/
    cout << "\n do nothing, in dsytrf() ";
} // end dsytrf

//--------------------------------------------------------------------
//void dsytri(char uplo,  int  n,  double  *a,  int  lda,  int
//             *ipivot, int *info)
void dsytri(char, int iter_used, matrix beta, int maxter, int *iwork,
            double *work, int &info)
{
    cout << "\n do nothing, in dsytri() ";
    // dsytri computes the inverse of a real  symmetric  indefinite
    // matrix  A  using  the  factorization  A  =  U*D*U**T  or A =
    // L*D*L**T computed by DSYTRF.

    // CHARACTER * 1 UPLO
    // INTEGER N, LDA, INFO
    // INTEGER IPIVOT(*)
    // DOUBLE PRECISION A(LDA,*), WORK(*)

    /**************************************************************
     UPLO (input)
               Specifies whether the details of the factorization
               are stored as an upper or lower triangular matrix.
               = 'U':  Upper triangular, form is A = U*D*U**T;
               = 'L':  Lower triangular, form is A = L*D*L**T.

     N (input) The order of the matrix A.  N >= 0.

     A (input/output)
               On entry, the block diagonal matrix D and the mul-
               tipliers  used to obtain the factor U or L as com-
               puted by DSYTRF.

               On exit, if INFO = 0, the (symmetric)  inverse  of
               the  original  matrix.   If  UPLO = 'U', the upper
               triangular part of the inverse is formed  and  the
               part of A below the diagonal is not referenced; if
               UPLO =  'L'  the  lower  triangular  part  of  the
               inverse  is  formed  and  the  part of A above the
               diagonal is not referenced.

     LDA (input)
               The leading dimension of  the  array  A.   LDA  >=
               max(1,N).

     IPIVOT (input)
               Details of the interchanges and the  block  struc-
               ture of D as determined by DSYTRF.

     WORK (workspace)
               dimension(N)

     INFO (output)
               = 0: successful exit
               < 0: if INFO = -i, the i-th argument had an  ille-
               gal value
               > 0: if INFO = i, D(i,i) = 0; the matrix is singu-
               lar and its inverse could not be computed.
    *****************************************************************/
} // end dsytri

//---------------------------------------------------------
double ddot(int n,
            double *x,
            int incx,
            double *y,
            int incy)
{
    // ddot compute the dot product of x and y where x  and  y  are
    // n-vectors.
    double prod;

    if (n < 1 || incx <= 0 || incy <= 0)
    {
        cerr << "\n error in ddot, n < 1 or incx <= 0 or incy <= 0, ";
        return 0;
    }

    prod = 0.0;

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        prod += x[ix] * y[iy];
    }

    return prod;
} // end ddot

//-----------------------------------------------------------------------
complex < double> ddot(int n,
                       complex < double> *x,
                       int incx,
                       complex < double> *y,
                       int incy)
{
    // zdotc compute the dot product of conjg(x) and y where x  and
    // y are n-vectors.
    complex < double> prod;

    if (n < 1 || incx <= 0 || incy <= 0)
    {
        cerr << "\n error in ddot, n < 1 or incx <= 0 or incy <= 0, ";
        return 0;
    }

    prod = complex < double>(0.0, 0.0);

    for (int ix = 0, iy = 0;ix < n && iy < n;ix += incx, iy += incy)
    {
        prod += conj(x[ix]) * y[iy];
    }

    return prod;
} // end zdotc

// -----------------------------------------------------------------------
complex < double> ddot(const ComplexMatrix &a,
                       int i,
                       complex < double> *y)
{
    //compute the dot product of i_th row of matrix a and y where
    // y are a vector.
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr)
    {
        cout << "\n error in ddot, nr or nc < 1 or i out of range ";
        return 0;
    }

    complex < double> z;

    z = complex < double> (0, 0);

    for (int k = 0; k < nc; k++)
    {
        z += conj(a(i, k)) * y[k];
    }

    return z;
}

//--------------------------------
double ddot(const matrix &a,
            int i,
            double  *y)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr)
    {
        cout << "\n error in ddot, nr or nc < 1 or i out of range ";
        return 0;
    }

    double z = 0;

    for (int k = 0; k < nc; k++)
    {
        z += a(i, k) * y[k];
    }

    return z;
}

double ddot(const matrix &a,
            int i,
            const matrix &b,
            int j)
{
    int nr, nc;
    nr = a.nr;
    nc = a.nc;

    if (nr <= 1 || nc <= 1 || i < 0 || i >= nr ||
            j < 0 || j >= nr)
    {
        cout << "\n error in ddot, nr or nc < 1 or i or j out of range ";
        return 0;
    }

    double z = 0;

    for (int k = 0; k < nc; k++)
    {
        z += a(i, k) * b(j, k);
    }

    return z;
}

double dnrm2(const int n, const double *x, const int incx)
{
    // compute Euclidean length (12 norm) of vector x,
    if (n < 0 || incx <= 0)
    {
        cerr << "\n error in dnrm2, n < 0 or incx <= 0, ";
        return 0;
    }
    if (n == 0)
    {
        return 0;
    }

    double norm2=0.0;
    for (int ix=0; ix<n; ix+=incx)
    {
        norm2 += x[ix] * x[ix];
    }

    return sqrt(norm2);
}

// i-row of matrix a
double dnrm2(const matrix &a,
             int i)
{
    int nc, nr;
    nc = a.nc;
    nr = a.nr;
    double norm2;
    norm2 = 0.0;

    if (nr < 1 || nc < 1 || i < 0 || i >= nr)
    {
        cerr << "\n error in dnrm2, nr or nc < 1 or i out of range ";
        return 0;
    }

    for (int j = 0;j < nc; j++)
    {
        norm2 += a(i, j) * a(i, j);
    }

    return sqrt(norm2);
}

void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           complex < double> alpha,
           const complex<double> *a,
           int lda,
           const ComplexMatrix &b,
           int ldb,
           complex < double> beta,
           complex<double> *c,
           int ldc)
{
	TITLE("myfunc5","zgemm1");
    //int nra=1;
	int nca=k;
	int nrb=b.nr;
	int ncb=b.nc;
	int nrc=1;
	int ncc=n;

    //double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        ij = 2*i;
        aux[ij] = a[i].real();
        aux[ij+1] = a[i].imag();
    }
    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * i;
            cux[ij]  = c[i].real();
            cux[ij+1] = c[i].imag();
        }
    }

	//WARNING_QUIT("reset zgemm","reset zgemm");
//    zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        ij = 2 *  i * ldc;
        c[i] = cux[ij] + complex< double>(0, 1) * cux[ij+1];
    }

    delete [] aux;
    delete [] bux;
    delete [] cux;
}


void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           complex < double> alpha,
           const complex<double> *a,
           int lda,
           const ComplexMatrix &b,
           int ldb,
           complex < double> beta,
           ComplexMatrix &c,
           int ldc)
{
    TITLE("myfunc5","zgemm2");
    //int nra = 1;
	int nca = k;
	int nrb = b.nr;
	int ncb = b.nc;
	int nrc = c.nr;
	int ncc = c.nc;

	//double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*nrc*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        ij = 2*i;
        aux[ij] = a[i].real();
        aux[ij+1] = a[i].imag();
    }
    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * nrc);
            cux[ij]  = c(j, i).real();
            cux[ij+1] = c(j, i).imag();
        }
    }

	//WARNING_QUIT("reset zgemm","reset zgemm");
    //zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2*(j+i*ldc);
            c(j,i) = cux[ij] + complex< double>(0,1) * cux[ij+1];
            //		cout<<setw(12)<<c(j,i);
        }
        //	cout<<endl;
    }

    delete [] aux;
    delete [] bux;
    delete [] cux;
}


// C = alpha * op(A) * op(B) + beta * C
void zgemm(char tra,
           char trb,
           int m,
           int n,
           int k,
           complex < double> alpha,
           const ComplexMatrix &a,
           int lda,
           const ComplexMatrix &b,
           int ldb,
           complex < double> beta,
           ComplexMatrix &c,
           int ldc)
{
//	TITLE("myfunc5","zgemm3");
    int nra, nca, nrb, ncb, nrc, ncc ;
    nra = a.nr,   nca = a.nc;
    nrb = b.nr,   ncb = b.nc;
    nrc = c.nr,   ncc = c.nc;

	//double alpha0[2] = {alpha.real(), alpha.imag()};
	//double beta0[2] = {beta.real(), beta.imag()};

    double *aux, *bux, *cux; //HLX: bug fixed on 12/20/2006
    aux = new double[2*nra*nca];
    bux = new double[2*nrb*ncb];
    cux = new double[2*nrc*ncc];

    int i;
    int j;
    int ij;

    for (i = 0; i < nca; i++)
    {
        for (j = 0; j < nra; j++)
        {
            ij = 2 * (j + i * nra);
            aux[ij] = a(j, i).real();
            aux[ij+1] = a(j, i).imag();
        }
    }

    for (i = 0; i < ncb; i++)
    {
        for (j = 0; j < nrb; j++)
        {
            ij = 2 * (j + i * nrb);
            bux[ij]  = b(j, i).real();
            bux[ij+1] = b(j, i).imag();
        }
    }


    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * nrc);
            cux[ij]  = c(j, i).real();
            cux[ij+1] = c(j, i).imag();
        }
    }

	//WARNING_QUIT("reset zgemm","reset zgemm");
    //zgemm_(&tra, &trb, &m, &n, &k, alpha0, aux, &lda, bux, &ldb, beta0, cux, &ldc);
	WARNING_QUIT("myfunc_5::zgemm","please don't ever use it again.");

    for (i = 0; i < ncc; i++)
    {
        for (j = 0; j < nrc; j++)
        {
            ij = 2 * (j + i * ldc);
            c(j, i) = cux[ij] + complex< double>(0, 1) * cux[ij+1];
        }
    }

    delete [] aux;

    delete [] bux;
    delete [] cux;
}//end zgemm


// C = alpha * op(A) * op(B) + beta * C
void dgemm(char tra, char trb, int m, int n, int k, double alpha,
           const matrix  a, int lda, const matrix b, int ldb, double beta,
           matrix &c, int ldc)
{
//    cout << "\n === ZGEMM() ===" << endl;
    matrix a1;

    if (tra == 'n' || tra == 'N')
        a1.create(a.nr, a.nc);
    else
        a1.create(a.nc, a.nr);

//    cout << "\n a1.nr = " << a1.nr
//	 << "   a1.nc = " << a1.nc << endl;

    c = beta * c;

    if (tra == 'n' || tra == 'N')
        a1 = alpha * a;
    else
        a1 = alpha * transpose(a);

    if (trb == 'n' || trb == 'N')
        c += a1 * b;
    else
        c += a1 * transpose(b);

//    a1.freemem();
//    cout << "\n end ZGEMM() " << endl;
    return;
}//end dgemm

int ILAENV(int ispec, char *name, char *opts,
           const int n1, const int n2, const int n3, const int n4)
{
    const int nb = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
    return nb;
}

void ZHPEV(int ,
           complex < double> *hp,
           double *e,
           ComplexMatrix &v,
           int ldh,
           int n,
           complex < double> *aux,
           int naux)
{
    cout << "\n do nothing, in ZHPEV() ";
}

complex < double> ZDOTU(int nstart,
                        complex < double>,
                        int	,
                        complex < double> *psi,
                        int npwx)
{
    cout << "\n do nothing in ZDOTU(), only return ZERO,";
    return ZERO;
}

void zgemv(char ,
           int ,
           int ,
           complex < double> alpha ,
           ComplexMatrix overlap,
           int ,
           complex < double> swfcatom ,
           int npwx,
           complex < double>  ,
           ComplexMatrix work, int)
{
    cout << "\n do nothing, in dgemv () ";
}

/*
void ZHEGVX(int itype,
            char jobz,
            char range ,
            char uplo ,
            const int n,
            const ComplexMatrix &a,
            const int lda,
            const ComplexMatrix &b,
            const int ldb,
            double vl,
            double vu,
            int il ,
            int iu,
            double abstol,
            int &m,
            double *w,
            ComplexMatrix &z,
            const int ldz,
            double *work,
            int lwork,
            double *rwork,
            int *iwork,
            int *ifail,
            int &info )
{
//	TITLE("myfunc5","ZHEGVX");
    double *aux, *bux, *zux;
    aux = new double[2*lda*n];//mohan fix + --> * 2007-10-22
    bux = new double[2*ldb*n];
    zux = new double[2*ldz*iu]; // mohan fix 2007-10-15
    int i, j;
	
	cout << "\n n = " << n;
	cout << "\n iu = " << iu << endl;

    for (i = 0;i < n;i++)
    {
        for (j = 0;j < lda;j++)
        {
            aux[2*(j+i*lda)]   = a(j, i).real();
            aux[2*(j+i*lda)+1] = a(j, i).imag();
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < ldb; j++)
        {
            bux[2*(j+i*ldb)]   = b(j, i).real();
            bux[2*(j+i*ldb)+1] = b(j, i).imag();
        }
    }

//	BLOCK_HERE("adf");
    zhegvx_(&itype, &jobz, &range, &uplo, &n, aux, &lda, bux, &ldb, &vl, &vu, &il, &iu, &abstol,
          &m, w, zux, &ldz, work, &lwork, rwork, iwork, ifail, &info);
	//BLOCK_HERE("adf2");

    for (i = 0;i < iu;i++)
    {
        for (j = 0;j < ldz;j++)
        {
            z(j, i) = zux[2*(j+i*ldz)] + complex< double>(0, 1) * zux[2*(j+i*ldz)+1];
        }
    }
    delete[] aux;
    delete[] bux;
    delete[] zux;
    return;
}
*/


