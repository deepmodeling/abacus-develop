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
 *  unit test of functions in cal_totn.cpp
 ***********************************************/

/**
 * - Tested functions in cal_totn.cpp:
 *   - cal_totn
 *     - calculate total charge density
 *   - induced_charge
 *     - calculate charges induced by the implicit solvent
 */

class cal_totn_test : public testing::Test
{
protected:
    
    surchem solvent_model;
};

TEST_F(cal_totn_test, cal_totn)
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
    complex<double>* TOTN = new complex<double>[npw];
    complex<double>* Porter_g = new complex<double>[npw];
    double *vloc = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(Porter_g,npw);

    for (int i=0; i<npw; i++)
    {
        Porter_g[i] = 0.1;
       
    }
    vloc[0] = -0.5593041647;
    vloc[1] = -0.3305673229;
    vloc[2] = -0.1228953775;
    for (int i=3; i<nrxx; i++)
    {
        vloc[i] = 0.1;     
    }

    solvent_model.cal_totn(GlobalC::ucell,GlobalC::rhopw,Porter_g,N,TOTN,vloc);

    EXPECT_NEAR(TOTN[0].real(),-0.0999496256,1e-10);
    EXPECT_NEAR(TOTN[0].imag(),-1.299621928166352e-7,1e-10);

    delete [] N;
    delete [] TOTN;
    delete [] Porter_g;
    delete [] vloc;
}

TEST_F(cal_totn_test, induced_charge)
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

    double fac;
    fac = ModuleBase::e2 * ModuleBase::FOUR_PI /(GlobalC::ucell.tpiba2 * GlobalC::rhopw->gg[0]);
    complex<double> delta_phi {-2.0347933860e-05,4.5900395826e-07};
    complex<double> induced_charge;
    induced_charge = -delta_phi/fac;

    EXPECT_NEAR(induced_charge.real(),6.2646421140e-05,1e-9);
    EXPECT_NEAR(induced_charge.imag(),-1.4135417993e-06,1e-9);
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
