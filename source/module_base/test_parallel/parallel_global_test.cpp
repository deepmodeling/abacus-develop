#ifdef __MPI
#include "mpi.h"
#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include <complex>
#include <string>
#include <cstring>

/************************************************
 *  unit test of functions in parallel_global.cpp
 ***********************************************/

/**
 * The tested functions are:
 *   i. Parallel_Global::split_diag_world(), which is
 *   used in David diagonalization in pw basis 
 *   calculation.
 *   ii. Parallel_Global::split_grid_world()
 *   iii. Parallel_Global::MyProd(std::complex<double> *in,std::complex<double> *inout,int *len,MPI_Datatype *dptr);
 *   iv. Parallel_Global::init_pools();
 *   v. Parallel_Global::divide_pools(void);
 */

TEST(ParaGlobal,SplitGrid)
{
	// NPROC is set to 4 in parallel_global_test.sh
	if(GlobalV::NPROC==4)
	{
		Parallel_Global::split_grid_world(2);
		EXPECT_EQ(GlobalV::GSIZE,2);
		if(GlobalV::MY_RANK==0) EXPECT_EQ(GlobalV::GRANK,0);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(GlobalV::GRANK,1);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(GlobalV::GRANK,0);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(GlobalV::GRANK,1);
		Parallel_Global::split_grid_world(4);
		EXPECT_EQ(GlobalV::GSIZE,1);
		if(GlobalV::MY_RANK==0) EXPECT_EQ(GlobalV::GRANK,0);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(GlobalV::GRANK,0);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(GlobalV::GRANK,0);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(GlobalV::GRANK,0);
	}
	else
	{
		Parallel_Global::split_grid_world(GlobalV::NPROC);
		EXPECT_EQ(GlobalV::GSIZE,1);
		EXPECT_EQ(GlobalV::GRANK,0);
	}
	//std::cout<<GlobalV::MY_RANK<<" "<<GlobalV::NPROC<<" ";
	//std::cout<<GlobalV::GRANK<<" "<<GlobalV::GSIZE<<std::endl;
}

TEST(ParaGlobal,SplitDiag)
{
	// NPROC is set to 4 in parallel_global_test.sh
	if(GlobalV::NPROC==4)
	{
		Parallel_Global::split_diag_world(2);
		EXPECT_EQ(GlobalV::DSIZE,2);
		if(GlobalV::MY_RANK==0) EXPECT_EQ(GlobalV::DRANK,0);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(GlobalV::DRANK,0);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(GlobalV::DRANK,1);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(GlobalV::DRANK,1);
		Parallel_Global::split_diag_world(4);
		EXPECT_EQ(GlobalV::DSIZE,4);
		if(GlobalV::MY_RANK==0) EXPECT_EQ(GlobalV::DRANK,0);
		if(GlobalV::MY_RANK==1) EXPECT_EQ(GlobalV::DRANK,1);
		if(GlobalV::MY_RANK==2) EXPECT_EQ(GlobalV::DRANK,2);
		if(GlobalV::MY_RANK==3) EXPECT_EQ(GlobalV::DRANK,3);
	}
	else
	{
		Parallel_Global::split_diag_world(GlobalV::NPROC);
		EXPECT_EQ(GlobalV::DSIZE,GlobalV::NPROC);
	}
	//std::cout<<GlobalV::MY_RANK<<" "<<GlobalV::NPROC<<" ";
	//std::cout<<GlobalV::DRANK<<" "<<GlobalV::DSIZE<<std::endl;
}

TEST(ParaGlobal, MyProd)
{
    std::complex<double> in[2] = {std::complex<double>(1.0, 2.0), std::complex<double>(-1, -2)};
    std::complex<double> inout[2] = {std::complex<double>(2.0, 1.0), std::complex<double>(-2, -1)};

    int len = 2;
    MPI_Datatype dptr = MPI_DOUBLE_COMPLEX;
    Parallel_Global::myProd(in, inout, &len, &dptr);
    EXPECT_EQ(inout[0], std::complex<double>(3.0, 3.0));
    EXPECT_EQ(inout[1], std::complex<double>(-3.0, -3.0));
}

TEST(ParaGlobal, InitPools)
{
    int nproc = 12;
    int kpar = 3;
    int stogroup = 3;
    int myrank = 5;
    GlobalV::NPROC = nproc;
    GlobalV::KPAR = kpar;
    GlobalV::NSTOGROUP = stogroup;
    GlobalV::MY_RANK = myrank;

    Parallel_Global::init_pools();
    EXPECT_EQ(GlobalV::NPROC_IN_STOGROUP, 4);
    EXPECT_EQ(GlobalV::MY_STOGROUP, 1);
    EXPECT_EQ(GlobalV::RANK_IN_STOGROUP, 1);
    EXPECT_EQ(GlobalV::MY_POOL, 0);
    EXPECT_EQ(GlobalV::RANK_IN_POOL, 1);
    EXPECT_EQ(GlobalV::NPROC_IN_POOL, 2);
    EXPECT_EQ(MPI_COMM_WORLD != STO_WORLD, true);
    EXPECT_EQ(STO_WORLD != POOL_WORLD, true);
    EXPECT_EQ(MPI_COMM_WORLD != PARAPW_WORLD, true);
}

TEST(ParaGlobal, DividePools)
{
    int nproc = 12;
    int kpar = 3;
    int stogroup = 3;
    int myrank = 5;
    GlobalV::NPROC = nproc;
    GlobalV::KPAR = kpar;
    GlobalV::NSTOGROUP = stogroup;
    GlobalV::MY_RANK = myrank;

    Parallel_Global::divide_pools();
    EXPECT_EQ(GlobalV::NPROC_IN_STOGROUP, 4);
    EXPECT_EQ(GlobalV::MY_STOGROUP, 1);
    EXPECT_EQ(GlobalV::RANK_IN_STOGROUP, 1);
    EXPECT_EQ(GlobalV::MY_POOL, 0);
    EXPECT_EQ(GlobalV::RANK_IN_POOL, 1);
    EXPECT_EQ(GlobalV::NPROC_IN_POOL, 2);
    EXPECT_EQ(MPI_COMM_WORLD != STO_WORLD, true);
    EXPECT_EQ(STO_WORLD != POOL_WORLD, true);
    EXPECT_EQ(MPI_COMM_WORLD != PARAPW_WORLD, true);
}

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif
