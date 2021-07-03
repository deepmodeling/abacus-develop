//==========================================================
// AUTHOR : fangwei, mohan
// DATE : 2009-11-08
//==========================================================
#include "parallel_global.h"
#include "parallel_common.h"
#include "parallel_reduce.h"
#include "../module_base/global_function.h"

#ifdef _OPENMP
#include <omp.h>					// Peize Lin add 2018-02-13
#endif

using namespace std;

#if defined __MPI
MPI_Datatype mpicomplex;
MPI_Op myOp;
MPI_Comm POOL_WORLD;
MPI_Comm PARAPW_WORLD; // qianrui add it for sto-dft 2021-4-14
MPI_Comm GRID_WORLD; // mohan add 2012-01-13z
MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

void Parallel_Global::myProd(complex<double> *in,complex<double> *inout,int *len,MPI_Datatype *dptr)
{
	for(int i=0;i<*len;i++)
	{
//		(*inout).real()=(*inout).real()+(*in).real();
//		(*inout).imag()=(*inout).imag()+(*in).imag();

		// mohan updat 2011-09-21
		(*inout)=complex<double>((*inout).real()+(*in).real(),(*inout).imag()+(*in).imag());

		in++;
		inout++;
	}
	return;
}
#endif


void Parallel_Global::split_diag_world(const int &diag_np)
{
#ifdef __MPI
	assert(diag_np>0);
	// number of processors in each 'grid group'.
	int* group_grid_np = new int[diag_np];
	ZEROS(group_grid_np, diag_np);
	// average processors in each 'grid group'
	int ave = NPROC/diag_np;
	// remain processors.
	int remain = NPROC - ave * diag_np;

	for(int i=0; i<diag_np; ++i)
	{
		group_grid_np[i] = ave;
		if(i<remain)
		{
			++group_grid_np[i];
		}
	}

	// color: same color will stay in same group.
	// key: rank in each fragment group.
	int color = -1;		// Peize Lin add initialization for compiler warning at 2020.01.31
	int key = -1;		// Peize Lin add initialization for compiler warning at 2020.01.31

	int np_now = 0;
	for(int i=0; i<diag_np; ++i)
	{
		np_now += group_grid_np[i];
		if(MY_RANK < np_now)
		{
			key = i;
			color = group_grid_np[i] - (np_now - MY_RANK);
			break;
		}
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, key, &DIAG_WORLD);
	MPI_Comm_rank(DIAG_WORLD, &DRANK);
	MPI_Comm_size(DIAG_WORLD, &DSIZE);
	DCOLOR=color;


	delete[] group_grid_np;
#else
	DCOLOR=0; //mohan fix bug 2012-02-04
	DRANK=0;
	DSIZE=1;
#endif
	return;
}



void Parallel_Global::split_grid_world(const int &diag_np)
{
#ifdef __MPI
	assert(diag_np>0);
	// number of processors in each 'grid group'.
	int* group_grid_np = new int[diag_np];
	ZEROS(group_grid_np, diag_np);
	// average processors in each 'grid group'
	int ave = NPROC/diag_np;
	// remain processors.
	int remain = NPROC - ave * diag_np;

	for(int i=0; i<diag_np; ++i)
	{
		group_grid_np[i] = ave;
		if(i<remain)
		{
			++group_grid_np[i];
		}
	}

	// color: same color will stay in same group.
	// key: rank in each fragment group.
	int color = -1;		// Peize Lin add initialization for compiler warning at 2020.01.31
	int key = -1;		// Peize Lin add initialization for compiler warning at 2020.01.31

	int np_now = 0;
	for(int i=0; i<diag_np; ++i)
	{
		np_now += group_grid_np[i];
		if(MY_RANK < np_now)
		{
			color = i;
			key = group_grid_np[i] - (np_now - MY_RANK);
			break;
		}
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, key, &GRID_WORLD);
	MPI_Comm_rank(GRID_WORLD, &GRANK);
	MPI_Comm_size(GRID_WORLD, &GSIZE);

	delete[] group_grid_np;
#else
	GRANK=0;  //mohan fix bug 2012-02-04
	GSIZE=1;
#endif
	return;
}

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

	#ifdef _OPENMP
//	omp_set_nested(true);					// Peize Lin add 2018-02-13
	#endif

//	MPI_Init(&argc,&argv);					// Peize Lin change 2018-07-12
	int provided;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
	if( provided != MPI_THREAD_FUNNELED )
		ofs_warning<<"MPI_Init_thread request "<<MPI_THREAD_FUNNELED<<" but provide "<<provided<<endl;
//----------------------------------------------------------
// int atoi ( const char * str );
// atoi : Convert string to int type
// atof : Convert string to double type
// atol : Convert string to long int type
//----------------------------------------------------------
//  NPOOL = atoi(argv[1]); // mohan abandon 2010-06-09

	// get the size --> NPROC
	// get the rank --> MY_RANK
	MPI_Comm_size(MPI_COMM_WORLD,&NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&MY_RANK);



	// for test
	for (int i=0; i<NPROC; i++)
    {
        if (MY_RANK == i)
        {
			if(i==0)
			{
				/*
				printf( "\n\e[33m%s\e[0m\n", " ===================================================");
				printf( "\e[33m%s\e[0m", "               WELCOME");
				printf( "\e[33m%s\e[0m", " TO");
				printf( "\e[33m%s\e[0m", " ESP");
				printf( "\e[33m%s\e[0m\n", " WORLD                 ");
				printf( "\e[33m%s\e[0m\n", " ===================================================");
				*/
				//xiaohui modify 2015-03-25
				/*
				cout << " *********************************************************" << endl;
				cout << " *                                                       *" << endl;
				cout << " *                  WELCOME TO MESIA                     *" << endl;
				cout << " *                                                       *" << endl;
				cout << " *       'Massive Electronic simulation based on         *" << endl;
				cout << " *        Systematically Improvable Atomic bases'        *" << endl;
				cout << " *                                                       *" << endl;
				cout << " *********************************************************" << endl;
				*/
				cout << " *********************************************************" << endl;
				cout << " *                                                       *" << endl;
				cout << " *                  WELCOME TO ABACUS                    *" << endl;
				cout << " *                                                       *" << endl;
				cout << " *            'Atomic-orbital Based Ab-initio            *" << endl;
				cout << " *                  Computation at UStc'                 *" << endl;
				cout << " *                                                       *" << endl;
                cout << " *          Website: http://abacus.ustc.edu.cn/          *" << endl;
                cout << " *                                                       *" << endl;
				cout << " *********************************************************" << endl;

				//cout << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    			time_t  time_now = time(NULL);
    			cout << " " << ctime(&time_now);
			}
//            cout << " PROCESSOR " << setw(4) << MY_RANK+1 << " IS READY." << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

	// This section can be chosen !!
	// mohan 2011-03-15
    if (MY_RANK != 0 )
    {
        //cout.rdbuf(NULL);
		cout.setstate(ios::failbit);//qianrui modify 2020-10-14
    }
	// end test

	MPI_Datatype block[2];
	block[0]=MPI_DOUBLE;
	block[1]=MPI_DOUBLE;

	int ac[2]={1,1};
	MPI_Aint dipc[2]={0,sizeof(double)};

	// MPI_Type_struct: create a struct datatype
	MPI_Type_create_struct(
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
