//==========================================================
// AUTHOR : fangwei, mohan
// DATE : 2009-11-08
//==========================================================
#include "parallel_global.h"

#ifdef __MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <thread>

#include "module_base/global_function.h"
#include "module_base/parallel_common.h"
#include "module_base/parallel_reduce.h"
#include "version.h"

#if defined __MPI
MPI_Comm POOL_WORLD;
MPI_Comm INTER_POOL;
MPI_Comm STO_WORLD;
MPI_Comm PARAPW_WORLD; // qianrui add it for sto-dft 2021-4-14
MPI_Comm GRID_WORLD; // mohan add 2012-01-13z
MPI_Comm DIAG_WORLD; // mohan add 2012-01-13

namespace Parallel_Global{
int mpi_number=0;
int omp_number=0;
}


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
	int provided = 0;
	MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
	if( provided != MPI_THREAD_MULTIPLE )
		GlobalV::ofs_warning<<"MPI_Init_thread request "<<MPI_THREAD_MULTIPLE<<" but provide "<<provided<<std::endl;
	// Peize Lin change 2022.08.08
	// MPI_THREAD_FUNNELED is enough for ABACUS. Using MPI_THREAD_SERIALIZED for elpa, using MPI_THREAD_MULTIPLE for libRI.
#else
	MPI_Init(&argc,&argv);					// Peize Lin change 2018-07-12
#endif //_OPENMP

    //  GlobalV::KPAR = atoi(argv[1]); // mohan abandon 2010-06-09

	// get the size --> GlobalV::NPROC
	// get the rank --> GlobalV::MY_RANK
	MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    // determining appropriate thread number for OpenMP
    const int max_thread_num = std::thread::hardware_concurrency(); // Consider Hyperthreading disabled.
#ifdef _OPENMP
    int current_thread_num = omp_get_max_threads();
#else
    int current_thread_num = 1;
#endif
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
    int process_num = 0, local_rank = 0;
    MPI_Comm_size(shmcomm, &process_num);
    MPI_Comm_rank(shmcomm, &local_rank);
    MPI_Comm_free(&shmcomm);
    mpi_number = process_num;
    omp_number = current_thread_num;
    if (current_thread_num * process_num > max_thread_num && local_rank==0)
    {
        std::stringstream mess;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%% WARNING: Total thread number(" << current_thread_num * process_num <<  ") " 
             << "is larger than hardware availability(" << max_thread_num << ")." << std::endl;
        mess << "%% WARNING: The results may be INCORRECT. Please be sure what you are doing." << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
        mess << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
		std::cerr << mess.str() << std::endl;
    }
    else if (current_thread_num * process_num < max_thread_num && local_rank==0)
    {
    	// only output info in local rank 0
        std::cerr << "WARNING: Total thread number on this node mismatches with hardware availability. "
            "This may cause poor performance."<< std::endl;
        std::cerr << "Info: Local MPI proc number: " << process_num << ","
                  << "OpenMP thread number: " << current_thread_num << ","
                  << "Total thread number: " << current_thread_num * process_num << ","
                  << "Local thread limit: " << max_thread_num << std::endl;
    }


    if (GlobalV::MY_RANK == 0)
    {
#ifdef VERSION
        const char* version = VERSION;
#else
        const char* version = "unknown";
#endif
#ifdef COMMIT_INFO
#include "commit.h"
        const char* commit = COMMIT;
#else
        const char* commit = "unknown";
#endif
        std::cout
            << "                                                                                     " << std::endl
            << "                              ABACUS " << version << std::endl
            << std::endl
            << "               Atomic-orbital Based Ab-initio Computation at UStc                    " << std::endl
            << std::endl
            << "                     Website: http://abacus.ustc.edu.cn/                             " << std::endl
            << "               Documentation: https://abacus.deepmodeling.com/                       " << std::endl
            << "                  Repository: https://github.com/abacusmodeling/abacus-develop       " << std::endl
            << "                              https://github.com/deepmodeling/abacus-develop         " << std::endl
            << "                      Commit: " << commit << std::endl
            << std::endl;
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
		std::cout.setstate(std::ios::failbit);//qianrui modify 2020-10-14
    }
	// end test
#endif //__MPI
    return;
}

#ifdef __MPI
void Parallel_Global::finalize_mpi()
{
	MPI_Comm_free(&POOL_WORLD);
    if (GlobalV::NPROC_IN_STOGROUP % GlobalV::KPAR == 0)
    {
        MPI_Comm_free(&INTER_POOL);
    }
	MPI_Comm_free(&STO_WORLD);
	MPI_Comm_free(&PARAPW_WORLD);
	MPI_Comm_free(&GRID_WORLD);
	MPI_Comm_free(&DIAG_WORLD);
	MPI_Finalize();
}
#endif

void Parallel_Global::init_pools(void)
{
#ifdef __MPI
//----------------------------------------------------------
// CALL Function : divide_pools
//----------------------------------------------------------
    Parallel_Global::divide_pools();

// for test
// turn on when you want to check the index of pools.
/*
    if (GlobalV::MY_RANK==0)
    {
        std::cout << "\n     " << std::setw(8) << "MY_RANK"
             << std::setw(8) << "MY_POOL"
             << std::setw(13) << "RANK_IN_POOL"
             << std::setw(6) << "NPROC"
             << std::setw(6) << "KPAR"
             << std::setw(14) << "NPROC_IN_POOL" << std::endl;
    }
    for (int i=0; i<GlobalV::NPROC; i++)
    {
        if (GlobalV::MY_RANK == i)
        {
            std::cout << " I'm " << std::setw(8) << GlobalV::MY_RANK
                 << std::setw(8) << GlobalV::MY_POOL
                 << std::setw(13) << GlobalV::RANK_IN_POOL
                 << std::setw(6) << GlobalV::NPROC
                 << std::setw(6) << GlobalV::KPAR
                 << std::setw(14) << GlobalV::NPROC_IN_POOL << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (GlobalV::MY_RANK != 0 )
    {
        std::cout.rdbuf(NULL);
    }
*/

    return;
#endif
}

#ifdef __MPI
void Parallel_Global::divide_pools(void)
{
    if (GlobalV::NPROC < GlobalV::KPAR)
    {
        std::cout<<"\n NPROC=" << GlobalV::NPROC << " KPAR=" << GlobalV::KPAR;
        std::cout<<"Error : Too many pools !"<<std::endl;
        exit(0);
    }

    // (1) per process in each stogroup
    if(GlobalV::NPROC%GlobalV::NSTOGROUP!=0)
    {
        std::cout<<"\n Error! NPROC="<<GlobalV::NPROC
        <<" must be divided evenly by BNDPAR="<<GlobalV::NSTOGROUP<<std::endl;
        exit(0);
    }
    GlobalV::NPROC_IN_STOGROUP = GlobalV::NPROC/GlobalV::NSTOGROUP;
    GlobalV::MY_STOGROUP = int(GlobalV::MY_RANK / GlobalV::NPROC_IN_STOGROUP);
    GlobalV::RANK_IN_STOGROUP = GlobalV::MY_RANK%GlobalV::NPROC_IN_STOGROUP;
    if (GlobalV::NPROC_IN_STOGROUP < GlobalV::KPAR)
    {
        std::cout<<"\n Error! NPROC_IN_BNDGROUP=" << GlobalV::NPROC_IN_STOGROUP 
            <<" is smaller than"<< " KPAR=" << GlobalV::KPAR<<std::endl;
        std::cout<<" Please reduce KPAR or reduce BNDPAR"<<std::endl;
        exit(0);
    }

    // (2) per process in each pool
    GlobalV::NPROC_IN_POOL = GlobalV::NPROC_IN_STOGROUP/GlobalV::KPAR;
    if (GlobalV::RANK_IN_STOGROUP < (GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR)*(GlobalV::NPROC_IN_POOL+1))
    {
        GlobalV::NPROC_IN_POOL++;
        GlobalV::MY_POOL = int(GlobalV::RANK_IN_STOGROUP / GlobalV::NPROC_IN_POOL);
        GlobalV::RANK_IN_POOL = GlobalV::RANK_IN_STOGROUP%GlobalV::NPROC_IN_POOL;
    }
    else
    {
        GlobalV::MY_POOL = int( (GlobalV::RANK_IN_STOGROUP-GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR) / GlobalV::NPROC_IN_POOL);
        GlobalV::RANK_IN_POOL = (GlobalV::RANK_IN_STOGROUP-GlobalV::NPROC_IN_STOGROUP%GlobalV::KPAR)%GlobalV::NPROC_IN_POOL;
    }
    



    int key = 1;
    MPI_Comm_split(MPI_COMM_WORLD,GlobalV::MY_STOGROUP,key,&STO_WORLD);

    //========================================================
    // MPI_Comm_Split: Creates new communicators based on
    // colors(2nd parameter) and keys(3rd parameter)
    // Note: The color must be non-negative or MPI_UNDEFINED.
    //========================================================
    MPI_Comm_split(STO_WORLD,GlobalV::MY_POOL,key,&POOL_WORLD);

    if (GlobalV::NPROC_IN_STOGROUP % GlobalV::KPAR == 0)
    {
        MPI_Comm_split(STO_WORLD, GlobalV::RANK_IN_POOL, key, &INTER_POOL);
    }

    int color = GlobalV::MY_RANK % GlobalV::NPROC_IN_STOGROUP;
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &PARAPW_WORLD);

    return;
}
#endif
