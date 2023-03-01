#ifndef __PWTEST
#define __PWTEST
#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
using namespace std;
extern int nproc_in_pool, rank_in_pool;
extern string precision_flag, device_flag;

class PWTEST: public testing::Test
{
public:
    static void SetUpTestCase()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"=    PW MODULE TEST START  ="<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
        }
    }
    static void TearDownTestCase()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"=     PW MODULE TEST END   ="<<"\033[0m"<<endl;
            cout<<"\033[32m"<<"============================"<<"\033[0m"<<endl;
        }
    }
    void SetUp()
    {
        cout<<"\033[32m"<<"[    CASE  ]"<<"\033[0m"<<" ";
    }
    void TearDown(){}
};

//memory.cpp depends on GlobalV::ofs_running and reduce_double_all
//GPU depends on reduce_double_pool
namespace GlobalV
{ 
    extern std::ofstream ofs_running;
}
#ifdef __MPI
#include "mpi.h"
extern MPI_Comm POOL_WORLD;
namespace Parallel_Reduce
{
    void reduce_double_all(double &object);
    void reduce_double_pool(double &object);
    void reduce_double_pool(float &object);
}
#endif
#endif