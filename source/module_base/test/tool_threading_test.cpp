#include "../tool_threading.h"
#include "gtest/gtest.h"
/************************************************
*  unit test of threading tool
***********************************************/

/**
* - Tested functions of class threading tool:
*   - TASK_DIST_1D:
*       - (template)Distributing 1d tasks by worker id (int and long long)
**/
class tool_threading_test : public testing::Test
{

};
TEST(tool_threading, TASK_DIST_1D_test_int)
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

TEST(tool_threading, TASK_DIST_1D_test_longlong)
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
