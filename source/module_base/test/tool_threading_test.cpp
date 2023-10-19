#include "../tool_threading.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <omp.h>
/************************************************
*  unit test of threading tool
***********************************************/

/**
* - Tested functions of class threading tool:
*   - TASK_DIST_1D:
*       - (template)Distributing 1d tasks by worker id (int and long long)
*   - BLOCK_TASK_DIST_1D:
*       - (template)Distributing 1d tasks by block_size and worker id (int and long long)
*   - OMP_PARALLE:
*       - Run functions in parallel mode
*   - TRY_OMP_PARALLEL:
*       - Run functions in parallel mode（Add the judgment statement to determine whether program is in parallel）
**/

// The meaning of the parameters used in the following tests
// nw: nworker
// iw: iworker
// nt: ntask
// st: start
// le: len
// bs: block_size

//Test function used in the following tests
void test_fun(int a,int b)
    {
        std::cout<<a<<std::endl;
        std::cout<<b;
        return ;
    }

TEST(ToolThreadingTEST, TastDist1DInt)
{
    int nw=1;
    int iw=1;
    int nt=16;
    int st;
    int le;
    ModuleBase::TASK_DIST_1D(nw,iw,nt,st,le);
    EXPECT_EQ(st,0);
    EXPECT_EQ(le,16);
    int nw2=4;
    ModuleBase::TASK_DIST_1D(nw2,iw,nt,st,le);
    EXPECT_EQ(st,4);
    EXPECT_EQ(le,4);
    int nw3=3;
    int nt3=17;
    ModuleBase::TASK_DIST_1D(nw3,iw,nt3,st,le);
    EXPECT_EQ(st,6);
    EXPECT_EQ(le,6);
}

TEST(ToolThreadingTEST, TastDist1DLonglong)
{
    long long nw=1;
    long long iw=1;
    long long nt=16;
    long long st;
    long long le;
    ModuleBase::TASK_DIST_1D(nw,iw,nt,st,le);
    EXPECT_EQ(st,0);
    EXPECT_EQ(le,16);
    long long nw2=4;
    ModuleBase::TASK_DIST_1D(nw2,iw,nt,st,le);
    EXPECT_EQ(st,4);
    EXPECT_EQ(le,4);
    long long nw3=3;
    long long nt3=17;
    ModuleBase::TASK_DIST_1D(nw3,iw,nt3,st,le);
    EXPECT_EQ(st,6);
    EXPECT_EQ(le,6);
}

TEST(ToolThreadingTEST, BlockTaskDist1DInt)
{
    int nw=1;
    int iw=1;
    int nt=16;
    int st;
    int le;
    int bs=1;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw,iw,nt,bs,st,le);
    EXPECT_EQ(st,0);
    EXPECT_EQ(le,16);
    int nw2=4;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw2,iw,nt,bs,st,le);
    EXPECT_EQ(st,4);
    EXPECT_EQ(le,4);
    int nw3=3;
    int nt3=17;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw3,iw,nt3,bs,st,le);
    EXPECT_EQ(st,6);
    EXPECT_EQ(le,6);
}

TEST(ToolThreadingTEST, BlockTaskDist1DLonglong)
{
    long long nw=1;
    long long iw=1;
    long long nt=16;
    long long st;
    long long le;
    long long bs=1;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw,iw,nt,bs,st,le);
    EXPECT_EQ(st,0);
    EXPECT_EQ(le,16);
    long long nw2=4;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw2,iw,nt,bs,st,le);
    EXPECT_EQ(st,4);
    EXPECT_EQ(le,4);
    long long nw3=3;
    long long nt3=17;
    ModuleBase:: BLOCK_TASK_DIST_1D(nw3,iw,nt3,bs,st,le);
    EXPECT_EQ(st,6);
    EXPECT_EQ(le,6);
}

TEST(ToolThreadingTEST, OmpParalle)
{
    std::string output1;
	testing::internal::CaptureStdout();
    ModuleBase::OMP_PARALLEL(test_fun);
	output1 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output1,testing::HasSubstr("1"));
    EXPECT_THAT(output1,testing::HasSubstr("0"));
}
#ifdef _OPENMP
TEST(ToolThreadingTEST, OmpParalleOpenmp)
{
    std::string output2;
	testing::internal::CaptureStdout();
    ModuleBase::OMP_PARALLEL(test_fun);
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(std::to_string(omp_get_num_threads())));
    EXPECT_THAT(output2,testing::HasSubstr(std::to_string(omp_get_thread_num())));
}
#undef _OPENMP
TEST(ToolThreadingTEST, TryOmpParalle)
{
    std::string output3;
	testing::internal::CaptureStdout();
    ModuleBase::TRY_OMP_PARALLEL(test_fun);
	output3 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output3,testing::HasSubstr("1"));
    EXPECT_THAT(output3,testing::HasSubstr("0"));
}
#define _OPENMP
TEST(ToolThreadingTEST, TryOmpParalleOpenmp)
{
    std::string output4;
	testing::internal::CaptureStdout();
    ModuleBase::TRY_OMP_PARALLEL(test_fun);
	output4= testing::internal::GetCapturedStdout();
	EXPECT_THAT(output4,testing::HasSubstr(std::to_string(omp_get_num_threads())));
    EXPECT_THAT(output4,testing::HasSubstr(std::to_string(omp_get_thread_num())));
}
#endif
