#ifdef __MPI
#include "../../src_parallel/parallel_global.h"
#include "mpi.h"
#include <iostream>
void setupmpi(int argc,char **argv,int &nproc, int &myrank)
{
    int provided;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
	if( provided != MPI_THREAD_FUNNELED )
		std::cout<<"MPI_Init_thread request "<<MPI_THREAD_FUNNELED<<" but provide "<<provided<<std::endl;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
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
}


void divide_pools(const int &nproc, const int &myrank, int &nproc_in_pool, int &kpar, int&mypool, int &rank_in_pool)
{  
    nproc_in_pool = nproc/kpar;
    if (myrank < (nproc%kpar)*(nproc_in_pool+1))
    {
        nproc_in_pool++;
    }

    int *nproc_pool = new int[kpar];
    int *startpro_pool = new int[kpar];
    for(int ip = 0 ; ip < kpar ; ++ip)
    {
        nproc_pool[ip] = 0;
        startpro_pool[ip] = 0;
    }

    for (int i=0; i<nproc; i++)
    {
        int j = i%kpar;
        nproc_pool[j]++;
    }

    // (3) To know start proc index in each pool.
    for (int i=1; i<kpar; i++)
    {
        startpro_pool[i]=startpro_pool[i-1]+nproc_pool[i-1];
    }

    // use 'myrank' to know 'mypool'.
    for (int i=0; i<kpar; i++)
    {
        if (myrank >= startpro_pool[i])
        {
            mypool=i;
        }
    }

    int key = 1;
    rank_in_pool = myrank-startpro_pool[mypool];

    MPI_Comm_split(MPI_COMM_WORLD,mypool,key,&POOL_WORLD);

    delete[] nproc_pool;
    delete[] startpro_pool;
    return;
}
void finishmpi()
{
    MPI_Comm_free(&POOL_WORLD);
    MPI_Type_free(&mpicomplex);
    MPI_Finalize();   
}
#endif