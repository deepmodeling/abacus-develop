//==========================================================
// AUTHOR : fangwei, mohan
// DATE : 2009-11-08
//==========================================================
#include "parallel_global.h"
#include "mpi.h"
#include "parallel_common.h"
#include "parallel_reduce.h"
#include "../module_base/global_function.h"
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#if defined __MPI
MPI_Datatype mpicomplex;
MPI_Op myOp;
MPI_Comm POOL_WORLD;
MPI_Comm STO_WORLD;
MPI_Comm PARAPW_WORLD; // qianrui add it for sto-dft 2021-4-14
MPI_Comm GRID_WORLD; // mohan add 2012-01-13z
MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

void Parallel_Global::myProd(std::complex<double> *in,std::complex<double> *inout,int *len,MPI_Datatype *dptr)
{
	for(int i=0;i<*len;i++)
	{
//		(*inout).real()=(*inout).real()+(*in).real();
//		(*inout).imag()=(*inout).imag()+(*in).imag();

		// mohan updat 2011-09-21
		(*inout)=std::complex<double>((*inout).real()+(*in).real(),(*inout).imag()+(*in).imag());

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
	ModuleBase::GlobalFunc::ZEROS(group_grid_np, diag_np);
	// average processors in each 'grid group'
	int ave = GlobalV::NPROC/diag_np;
	// remain processors.
	int remain = GlobalV::NPROC - ave * diag_np;

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
		if(GlobalV::MY_RANK < np_now)
		{
			key = i;
			color = group_grid_np[i] - (np_now - GlobalV::MY_RANK);
			break;
		}
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, key, &DIAG_WORLD);
	MPI_Comm_rank(DIAG_WORLD, &GlobalV::DRANK);
	MPI_Comm_size(DIAG_WORLD, &GlobalV::DSIZE);
	GlobalV::DCOLOR=color;


	delete[] group_grid_np;
#else
	GlobalV::DCOLOR=0; //mohan fix bug 2012-02-04
	GlobalV::DRANK=0;
	GlobalV::DSIZE=1;
#endif
	return;
}



void Parallel_Global::split_grid_world(const int &diag_np)
{
#ifdef __MPI
	assert(diag_np>0);
	// number of processors in each 'grid group'.
	int* group_grid_np = new int[diag_np];
	ModuleBase::GlobalFunc::ZEROS(group_grid_np, diag_np);
	// average processors in each 'grid group'
	int ave = GlobalV::NPROC/diag_np;
	// remain processors.
	int remain = GlobalV::NPROC - ave * diag_np;

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
		if(GlobalV::MY_RANK < np_now)
		{
			color = i;
			key = group_grid_np[i] - (np_now - GlobalV::MY_RANK);
			break;
		}
	}

	MPI_Comm_split(MPI_COMM_WORLD, color, key, &GRID_WORLD);
	MPI_Comm_rank(GRID_WORLD, &GlobalV::GRANK);
	MPI_Comm_size(GRID_WORLD, &GlobalV::GSIZE);

	delete[] group_grid_np;
#else
	GlobalV::GRANK=0;  //mohan fix bug 2012-02-04
	GlobalV::GSIZE=1;
#endif
	return;
}

void Parallel_Global::read_mpi_parameters(int argc,char **argv)
{
#ifdef __MPI
#ifdef _OPENMP
	int provided;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_SERIALIZED,&provided);
	if( provided != MPI_THREAD_SERIALIZED )
		GlobalV::ofs_warning<<"MPI_Init_thread request "<<MPI_THREAD_SERIALIZED<<" but provide "<<provided<<std::endl;
	// Peize Lin change 2021.11.24
	// MPI_THREAD_FUNNELED is enough for ABACUS. Using MPI_THREAD_SERIALIZED for elpa.
#else
	MPI_Init(&argc,&argv);					// Peize Lin change 2018-07-12
#endif //_OPENMP

    //  GlobalV::KPAR = atoi(argv[1]); // mohan abandon 2010-06-09

	// get the size --> GlobalV::NPROC
	// get the rank --> GlobalV::MY_RANK
	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    // determining appropriate thread number for OpenMP
#ifdef _OPENMP
    const int max_thread_num = std::thread::hardware_concurrency(); // Consider Hyperthreading disabled.
    int current_thread_num = omp_get_max_threads();
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    int process_num, local_rank;
    MPI_Comm_size(shmcomm, &process_num);
    MPI_Comm_rank(shmcomm, &local_rank);
    MPI_Comm_free(&shmcomm);
    int desired_thread_num = max_thread_num / process_num;
    if (desired_thread_num != current_thread_num && current_thread_num == max_thread_num)
    {
        // OpenMP thread num not set
        omp_set_num_threads(desired_thread_num);
        current_thread_num = omp_get_max_threads();
    }
    if (current_thread_num * process_num != max_thread_num && local_rank==0)
    {
		// only output info in local rank 0
        std::cerr << "WARNING: Total thread number on this node mismatches with hardware availability. "
            "This may cause poor performance."<< std::endl;
        std::cerr << "Info: Local MPI proc number: " << process_num << ","
                  << "OpenMP thread number: " << current_thread_num << ","
                  << "Total thread number: " << current_thread_num * process_num << ","
                  << "Local thread limit: " << max_thread_num << std::endl;
    }
#endif

    if (GlobalV::MY_RANK == 0)
    {
        std::cout << " *********************************************************" << std::endl;
        std::cout << " *                                                       *" << std::endl;
        std::cout << " *                  WELCOME TO ABACUS v3.0               *" << std::endl;
        std::cout << " *                                                       *" << std::endl;
        std::cout << " *            'Atomic-orbital Based Ab-initio            *" << std::endl;
        std::cout << " *                  Computation at UStc'                 *" << std::endl;
        std::cout << " *                                                       *" << std::endl;
        std::cout << " *          Website: http://abacus.ustc.edu.cn/          *" << std::endl;
        std::cout << " *                                                       *" << std::endl;
        std::cout << " *********************************************************" << std::endl;
        time_t time_now = time(NULL);
        std::cout << " " << ctime(&time_now);
    }

    // for test
    /*
	for (int i=0; i<GlobalV::NPROC; i++)
    {
        if (GlobalV::MY_RANK == i)
        {
            std::cout << " PROCESSOR " << std::setw(4) << GlobalV::MY_RANK+1 << " IS READY." << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    */

	// This section can be chosen !!
	// mohan 2011-03-15
    if (GlobalV::MY_RANK != 0 )
    {
        //std::cout.rdbuf(NULL);
		std::cout.setstate(ios::failbit);//qianrui modify 2020-10-14
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

#endif //__MPI
    return;
}

#ifdef __MPI
void Parallel_Global::finalize_mpi()
{
	MPI_Comm_free(&POOL_WORLD);
	MPI_Comm_free(&STO_WORLD);
	MPI_Comm_free(&PARAPW_WORLD);
	MPI_Comm_free(&GRID_WORLD);
	MPI_Comm_free(&DIAG_WORLD);
	MPI_Finalize();
}
#endif
