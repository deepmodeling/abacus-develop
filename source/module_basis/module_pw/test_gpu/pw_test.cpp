#ifdef __MPI
#include "mpi.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fftw3.h"
#include "pw_test.h"
using namespace std;

int nproc_in_pool, rank_in_pool;
string precision_flag, device_flag;

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
    precision_flag = "double";
    device_flag = "cpu";
    nproc_in_pool = kpar = 1;
    rank_in_pool = 0;
#ifdef __MPI
    int nproc;
    int myrank;
    int mypool;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#endif
#ifdef _OPENMP
    // ref: https://www.fftw.org/fftw3_doc/Usage-of-Multi_002dthreaded-FFTW.html
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    int result = 0;
    testing::AddGlobalTestEnvironment(new TestEnv);
    testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
#ifdef __MPI
   MPI_Finalize();
#endif
#ifdef _OPENMP
	fftw_cleanup_threads();
#endif
    return result;
}
