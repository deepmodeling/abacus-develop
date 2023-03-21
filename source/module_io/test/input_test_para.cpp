#include "gtest/gtest.h"
#include "module_base/global_variable.h"
#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of Input::bcast
 ***********************************************/

/**
 * - Tested Functions:
 *   - bcast()
 *     - bcast input parameters to all processes
 */

#define private public
#include "module_io/input.h"

class InputParaTest : public ::testing::Test
{
protected:
};

#ifdef __MPI
TEST_F(InputParaTest,Bcast)
{
	if(GlobalV::MY_RANK==0)
	{
		INPUT.Default();
	}
	INPUT.Bcast();
	if(GlobalV::MY_RANK!=0)
	{
		EXPECT_EQ(INPUT.suffix,"ABACUS");
	}
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
#undef private
