//column-circle decomposition
#include "pzgseps.h"
#include "src_parallel/parallel_reduce.h"
#include "module_base/lapack_connector.h"

#include "pzheg2st.h"
#include "pdstebz.h"
#include "pzhetrd.h"
#include "pzsteiz.h"
#include "pzt2s.h"
#include "pzst2g.h"

void pzgseps(MPI_Comm comm_2D,int n,int nb,int &egnum, complex<double> *A, complex<double> *B,complex<double> *Z,double *eigen,
             LocalMatrix LM, char uplo, int &loc_size,int &loc_pos)
{
    TITLE("Parrallel_Diago","pzgseps");
    int i,j,count=0,k=0,incx=1;
    int color,key,err;

    int dim[2],period[2],coord[2],NPROCS,IAM;
    complex<double> norm[n];
    double tol=0.00000001;

    MPI_Cart_get(comm_2D,2,dim,period,coord);
    IAM=coord[0]*dim[1]+coord[1];
    NPROCS=dim[0]*dim[1];

    double* diag = new double[n];
    double* off_diag = new double[n];
    double* D = new double[n];
    double* E = new double[n];
    complex<double>* Z1 = new complex<double>[loc_size*n];

	int wrong_input = 0;

	/*
	ofs_running << "\n nb = " << nb;
	ofs_running << "\n n = " << n;
	ofs_running << "\n LocalMatrix A information : ";
	ofs_running << "\n col_num = " << LM.col_num;
	ofs_running << "\n row_num = " << LM.row_num;
	ofs_running << "\n col_pos = " << LM.col_pos;
	ofs_running << "\n row_pos = " << LM.row_pos;
	ofs_running << "\n row_b   = " << LM.row_b;
	ofs_running << "\n col_b   = " << LM.col_b;
	*/

	/*
	for(int i=0; i<LM.row_num; i++)
	{
		ofs_running << "\n row_set[" << i << "] = " << LM.row_set[i];
	}

	for(int i=0; i<LM.col_num; i++)
	{
		ofs_running << "\n col_set[" << i << "] = " << LM.col_set[i];
	} 
	 */
//	ofs_running << "\n LocalMatrix Information done." << endl;

	if(LM.col_num==0 || LM.row_num==0)
	{
		wrong_input = 1;
	}

	Parallel_Reduce::reduce_int_all( wrong_input );

	//ofs_running << "\n wrong_input = " << wrong_input << endl;
	if(wrong_input > 0)
	{
		ofs_running << "\n col_num == 0 || row_num == 0" << endl;
		WARNING_QUIT("pzgseps","col_num == 0 || row_num == 0");
	}

	pzheg2st(uplo,comm_2D,A,LM,B,LM,nb,n);
	
//	for(int yang=0;yang<LM.col_num*LM.row_num;yang++)
//	{
//		printf("A[%d]=%lf,%lf ",yang,A[yang].real(),A[yang].imag());
//	}
	

//	for(int yang=0;yang<n;yang++)
//	{
//		printf("Z[%d]=%lf+i%lf,norm[%d]=%lf+i%lf ",yang,Z[yang].real(),Z[yang].imag(),yang,norm[yang].real(),norm[yang].imag());
//	}

	pzhetrd(comm_2D,LM,n,nb,A,diag,off_diag,norm,uplo);

//	for (i=0;i<n;i++)
//		printf("\neigen[%d]=%lf,diag[%d]=%lf,off_diag[%d]=%lf ",i,eigen[i],i,diag[i],i,off_diag[i]);

	LapackConnector::copy(n,diag,incx,D,incx);
	LapackConnector::copy(n,off_diag,incx,E,incx);

	pdstebz(comm_2D,diag,off_diag,eigen,n);


	loc_size=n/NPROCS;

	if (n%NPROCS==0)
	{
		loc_pos=loc_size*IAM;
	}
	else
	{
		if (IAM<n%NPROCS)
		{
			// mohan update 2010-12-29
			loc_pos=(loc_size+1)*IAM;
			//loc_pos=loc_size*IAM;
		}
		else
		{
			loc_pos=(loc_size+1)*(n%NPROCS)+loc_size*(IAM-n%NPROCS);
		}
	}
	if (IAM<n%NPROCS)
	{
		++loc_size;
	}

//	ofs_running << "\n loc_size = " << loc_size << endl;

//	printf("loc_pos=%d\n",loc_pos);
	pzsteiz(n,D,E,&eigen[loc_pos],Z,loc_size);
//	printf("loc_pos=%d\n",loc_pos);

//	printf("n=%d nb=%d loc_size=%d uplo=%c\n",n,nb,loc_size,uplo);

	pzt2s(comm_2D,n,nb,A,Z,LM,loc_size,norm,uplo);
	if (uplo=='l'||(uplo=='L'))
	{
		for(i=0;i<LM.row_num;i++)
			for(j=0;j<LM.col_num;j++)
				B[i*LM.col_num+j]=conj(B[i*LM.col_num+j]);
		pzst2g(comm_2D,nb,n,B,Z,LM,loc_size,uplo);
	}

//	ofs_running << "\n local_size = " << loc_size << endl;
	if (uplo=='u'||(uplo=='U'))
	{
		for (i=0; i<loc_size; i++)
		{
			for (j=0; j<n; j++)
			{
				Z1[j*loc_size+i]=Z[i*n+j];
//				cout << " i=" << i << " j=" << j << " Z=" << Z[i*n+j] << endl;
			}
		}
		pzst2g(comm_2D,nb,n,B,Z1,LM,loc_size,uplo);
		for (i=0; i<loc_size; i++)
		{
			for (j=0; j<n; j++)
			{
				Z[i*n+j]=Z1[i*n+j];
			}
		}
	}
	egnum=n;

	delete[] diag;
	delete[] off_diag;
	delete[] D;
	delete[] E;
	delete[] Z1;

//	ofs_running << "\n Finish." << endl;
	return;
}
