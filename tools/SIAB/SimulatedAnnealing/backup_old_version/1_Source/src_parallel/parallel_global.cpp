//==========================================================
// AUTHOR : fangwei, mohan
// DATE : 2009-11-08
//==========================================================
#include "parallel_global.h"
#include "parallel_common.h"
#include "parallel_reduce.h"
#include "../src_spillage/tools.h"

using namespace std;

#if defined __MPI
MPI_Datatype mpicomplex;
MPI_Op myOp;
MPI_Comm POOL_WORLD;

void Parallel_Global::myProd(complex<double> *in,complex<double> *inout,int *len,MPI_Datatype *dptr)
{
	for(int i=0;i<*len;i++)
	{
		(*inout).real()=(*inout).real()+(*in).real();
		(*inout).imag()=(*inout).imag()+(*in).imag();
		in++;
		inout++;
	}	
	return;
}
#endif

void Parallel_Global::read_mpi_parameters(int argc,char **argv)
{
#if defined __MPI
//for test
/*
    cout << "\n Hello! Test MPI NOW : argc = "<<argc<<endl;
    for(int i=0;i<argc;i++)
    {
        cout<<"\n argv["<<i<<"]="<<argv[i];
    }
    cout<<endl;
*/
    MPI_Init(&argc,&argv);
//----------------------------------------------------------
// int atoi ( const char * str );
// atoi : Convert string to int type
// atof : Convert string to double type
// atol : Convert string to long int type
//----------------------------------------------------------
//  KPAR = atoi(argv[1]); // mohan abandon 2010-06-09

	// get the size --> NPROC
	// get the rank --> MY_RANK
	MPI_Comm_size(MPI_COMM_WORLD,&NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&MY_RANK);

	for(int i=0; i<NPROC; i++)
	{
		if(i==MY_RANK)
		{
			cout << " My rank is " << MY_RANK << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	MPI_Datatype block[2];
	block[0]=MPI_DOUBLE;
	block[1]=MPI_DOUBLE;
	
	int ac[2]={1,1};
	MPI_Aint dipc[2]={0,sizeof(double)};

	// MPI_Type_struct: create a struct datatype
	MPI_Type_struct(
	2,// count: number of blocks(integer)
	ac,//number of element in each block(array)
	dipc,//byte displacement of each block
	block,//type of element in each block(array of handles to datatype objects)
	&mpicomplex);//new type
	
	MPI_Type_commit(&mpicomplex);
	MPI_Op_create((MPI_User_function *)Parallel_Global::myProd,1,&myOp);
	
#endif
    return;
}
