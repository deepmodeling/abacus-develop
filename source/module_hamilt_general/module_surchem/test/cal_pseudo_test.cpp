#ifdef __MPI
#include "mpi.h"
#include "../../../module_base/parallel_global.h"
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>
#include <fstream>
#include "../surchem.h"
#include "setcell.h"

/************************************************
 *  unit test of functions in cal_pseudo.cpp
 ***********************************************/

/**
 * - Tested functions in cal_pseudo.cpp:
 *   - gauss_charge
 *     - calculate the gaussian pseudo core charge
 *   - cal_pseudo
 *     - add pseudo charge to rho
**/

class cal_pseudo_test : public testing::Test
{
protected:
    
    surchem solvent_model;

};


TEST_F(cal_pseudo_test, gauss_charge)
{   
    std::string precision_flag, device_flag;
    precision_flag = "double";
    device_flag = "cpu";

    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    GlobalC::rhopw = &pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    bool gamma_only;

    Setcell::setupcell(GlobalC::ucell);

    wfcecut = 80;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;

    //init
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK); 
    Parallel_Global::split_diag_world(GlobalV::NPROC);
    Parallel_Global::split_grid_world(GlobalV::NPROC);
    MPI_Comm_split(MPI_COMM_WORLD,0,1,&POOL_WORLD); //in LCAO kpar=1
#endif

#ifdef __MPI
    GlobalC::rhopw ->initmpi(1, 0, POOL_WORLD);
#endif
    //pwtest.initgrids(lat0,latvec,wfcecut);

    GlobalC::rhopw ->initgrids(GlobalC::ucell.lat0,GlobalC::ucell.latvec,wfcecut);

    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();
    GlobalC::rhopw ->collect_uniqgg();

    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;
    complex<double> *N = new complex<double>[npw];
    ModuleBase::GlobalFunc::ZEROS(N, npw);

    Structure_Factor sf;
    sf.nbspline = -1;

    solvent_model.gauss_charge(GlobalC::ucell,GlobalC::rhopw,N,&sf);

    EXPECT_NEAR(N[14].real(),0.002,1e-9);
    EXPECT_NEAR(N[16].real(),-0.001573534,1e-9);

    delete [] N;

}


TEST_F(cal_pseudo_test, cal_pseudo)
{   
    std::string precision_flag, device_flag;
    precision_flag = "double";
    device_flag = "cpu";

    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    GlobalC::rhopw = &pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    bool gamma_only;

    //--------------------------------------------------
    wfcecut = 80;
    gamma_only = false;
    int distribution_type = 1;
    bool xprime = false;
    //--------------------------------------------------

    //init
#ifdef __MPI
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK); 
    Parallel_Global::split_diag_world(GlobalV::NPROC);
    Parallel_Global::split_grid_world(GlobalV::NPROC);
    MPI_Comm_split(MPI_COMM_WORLD,0,1,&POOL_WORLD); //in LCAO kpar=1
#endif

#ifdef __MPI
    GlobalC::rhopw ->initmpi(1, 0, POOL_WORLD);
#endif
    //pwtest.initgrids(lat0,latvec,wfcecut);

    GlobalC::rhopw ->initgrids(GlobalC::ucell.lat0,GlobalC::ucell.latvec,wfcecut);
    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();
    GlobalC::rhopw ->collect_uniqgg();

    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;
    
    Structure_Factor sf;
    sf.nbspline = -1;

    complex<double>* Porter_g = new complex<double>[npw];
    ModuleBase::GlobalFunc::ZEROS(Porter_g,npw);
    for (int i=0; i<npw; i++)
    {
       Porter_g[i] = 0.1;
    }
    
    complex<double>* PS_TOTN = new complex<double>[npw];
    solvent_model.cal_pseudo(GlobalC::ucell,GlobalC::rhopw,Porter_g,PS_TOTN,&sf);

    EXPECT_NEAR(PS_TOTN[16].real(),0.098426466,1e-9);
    EXPECT_NEAR(PS_TOTN[14].real(),0.102,1e-9);

    delete [] Porter_g;
    delete [] PS_TOTN;
}


int main(int argc, char **argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}


