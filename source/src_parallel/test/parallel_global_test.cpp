#ifdef __MPI
#include "mpi.h"
#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#include "src_parallel/parallel_global.h"
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
