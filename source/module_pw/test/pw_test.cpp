#ifdef __MPI
#include "test_tool.h"
#include "mpi.h"
#endif
#include "../../module_base/global_variable.h"
#include "fftw3.h"
#include "pw_test.h"
using namespace std;

int nproc_in_pool, rank_in_pool;

class TestEnv : public testing::Environment 
{
public:
    virtual void SetUp()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"[ SET UP TESTS ]"<<"\033[0m"<<endl;
            cout<<"\033[32m[ "<<nproc_in_pool<<" processors are used."<<" ]\033[0m"<<endl;
        }
    }
    virtual void TearDown()
    {
        if(rank_in_pool == 0)
        {
            cout<<"\033[32m"<<"[ TEAR DOWN TESTS ]"<<"\033[0m"<<endl;
        }
    }
};



int main(int argc, char **argv) 
{
    
    int kpar;
    kpar = 1;
#ifdef __MIX_PRECISION
    //Temporary, pw_basis should not contains global variables
    GlobalV::precision_flag = "single";
#endif
#ifdef __MPI
    int nproc, myrank,mypool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#else
    nproc_in_pool = kpar = 1;
    rank_in_pool = 0;
#endif
    int result = 0;
    testing::AddGlobalTestEnvironment(new TestEnv);
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif  
    return result;
}
