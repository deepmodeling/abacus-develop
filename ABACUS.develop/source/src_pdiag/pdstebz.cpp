#include "pdstebz.h"
#include "module_base/lapack_connector.h"

void pdstebz(MPI_Comm comm ,double *D,double *E,double *eigen,int N)
/*
 *  PSEPS routine (version 1.0) --
 *  Computer Network Information Center, CAS.
 *  December 15, 2004
 *
 *  Purpose
 *  ========
 *  PSSTEBZ computes the eigenvalues of a symmetric tridiagonal matrix in
 *  parallel.The subroutine PSTEBZ computes the specified eigenvalues of T
 *  by divide-conquer algorithm based on rank-1 modification + laguerre iteration
 *
 *  Arguments
 *  =========
 *  N          (input) integer
 *             The order of the matrix A N>=0
 *  comm       (input)mpi communicator
 *  myid       (local input)rank number of this processor
 *  size       (input)total number of processor
 *  D          (local input) double precision array, dimension(n)
 *             The diagonal elements of the tridiagonal symmetric matrix
 *  E          (local input) double precision array, dimension(n-1)
 *             The off-diagonal elements of the tridiagonal matrix
 *  loc_eig    (output) double precision  array, dimension(n)
 *             Stores the eigenvalue required
 */
{
	TITLE("Parallel_Diago","pdstebz");
    /*Array A contains main diagonal element;Array B contains offdiagonal element*/
    int size,myid;

    int incx=1,incy=1;
    int dest,pmid,pstart,pend,info=0;
    MPI_Status status;
    double t1,t2;

    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&myid);

    /*Call dsterf in LAPACK ,where dsterf can solve eig by QR algorithm*/
    dsterf_(&N,D,E,&info);

	/*
	for(int i=0; i<N; i++)
	{
		eigen[i] = D[i];
		cout << " D[" << i << "]=" << D[i] << endl;
	}
	*/
    
	LapackConnector::copy(N,D,incx,eigen,incy);
	return;
}
