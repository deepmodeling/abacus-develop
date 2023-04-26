#include "../pw_basis.h"
#include "test_tool.h"
#include "module_base/parallel_global.h"
#include "mpi.h"
#include "module_base/global_function.h"
#include "module_base/constants.h"
#include "pw_test.h"
using namespace std;

void create_pools(const int totnproc, const int myrank, const int nproc)
{  
    int mypool = 1;
    if(myrank < nproc) mypool = 0;
    MPI_Comm_split(MPI_COMM_WORLD,mypool,1,&POOL_WORLD);
    return;
}

int main(int argc, char **argv)
{
    int totnproc, myrank;
    setupmpi(argc,argv,totnproc, myrank);
    const int nn = totnproc;
    int npw_per_ref[nn][nn];
    int nst_per_ref[nn][nn];
    int totnpw_ref;
    int totnst_ref;
    int nx_ref;
    int fftnx_ref;
    int ny_ref;
    int fftny_ref;
    int nz_ref;
    if(myrank == 0) cout<<"Generating ref..."<<endl;
    //--------------------------------------------------
    ModuleBase::Matrix3 latvec(1,2,0,2,1,1,0,0,3);
    bool gamma_only = true;
    double wfcecut = 70;
    double lat0 = 5;
    int distribution_type = 2;
    bool xprime = true;
    //--------------------------------------------------
    for(int nproc = 1 ; nproc <= totnproc ; ++nproc)
    {
      create_pools(totnproc, myrank, nproc);
      if(myrank < nproc)  
      {
        ModulePW::PW_Basis pwtest;
#ifdef __MPI
        pwtest.initmpi(nproc, myrank, POOL_WORLD);
#endif
        pwtest.initgrids(lat0, latvec, wfcecut);
        pwtest.initparameters(gamma_only, wfcecut, distribution_type,xprime);
        pwtest.setuptransform();
        pwtest.collect_local_pw();
        pwtest.collect_uniqgg();
        ModuleBase::Matrix3 GT,G,GGT;
        GT = latvec.Inverse();
	    G  = GT.Transpose();
	    GGT = G * GT;
        double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / lat0 / lat0;
        double ggecut = wfcecut / tpiba2;


        int *npw_per = nullptr;
        if(myrank == 0)
        {
            npw_per = new int [nproc];
        }

        MPI_Gather(&pwtest.npw,1,MPI_INT,npw_per,1,MPI_INT,0,POOL_WORLD);

        if(myrank == 0)
        {
            for(int ip = 0 ; ip < nproc ; ++ip)
            {
                npw_per_ref[nproc-1][ip] = npw_per[ip];
                nst_per_ref[nproc-1][ip] = pwtest.nst_per[ip];
            }
            delete []npw_per;
        }

        totnpw_ref = pwtest.npwtot;
        totnst_ref = pwtest.nstot;
        nx_ref = pwtest.nx;
        fftnx_ref = pwtest.fftnx;
        ny_ref = pwtest.ny;
        fftny_ref=pwtest.fftny;
        nz_ref =pwtest.nz;
      }
        MPI_Comm_free(&POOL_WORLD);
    }
    if(myrank == 0)
    {
        cout<<"const int totnpw_ref = "<<totnpw_ref<<";"<<endl;
        cout<<"const int totnst_ref = "<<totnst_ref<<";"<<endl;
        cout<<"const int nx_ref = "<<nx_ref<<";"<<endl;
        cout<<"const int fftnx_ref = "<<fftnx_ref<<";"<<endl;
        cout<<"const int ny_ref = "<<ny_ref<<";"<<endl;
        cout<<"const int fftny_ref = "<<fftny_ref<<";"<<endl;
        cout<<"const int nz_ref = "<<nz_ref<<";"<<endl;


        cout<<endl;

        cout<<"int npw_per_ref["<<totnproc<<"]["<<totnproc<<"]={"<<endl;
        for(int ip = 0 ; ip < totnproc ; ++ip)
        {
            cout<<"    {";
            for(int j = 0; j < ip+1 ; ++j)
            {
                cout<<npw_per_ref[ip][j];
                if(j<ip) cout<<",";
            }
            if(ip < totnproc - 1)   cout<<"},"<<endl;
            else                    cout<<"}"<<endl;
        }
        cout<<"};"<<endl;

        cout<<endl;

        cout<<"int nst_per_ref["<<totnproc<<"]["<<totnproc<<"]={"<<endl;
        for(int ip = 0 ; ip < totnproc ; ++ip)
        {
            cout<<"    {";
            for(int j = 0; j < ip+1 ; ++j)
            {
                cout<<nst_per_ref[ip][j];
                if(j<ip) cout<<",";
            }
            if(ip < totnproc - 1)   cout<<"},"<<endl;
            else                    cout<<"}"<<endl;
        }
        cout<<"};"<<endl;
        
    }

    MPI_Finalize(); 


    return 0;  
}