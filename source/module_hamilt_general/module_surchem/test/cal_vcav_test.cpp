#ifdef __MPI
#include "mpi.h"
#include "../../../module_base/parallel_global.h"
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>
#include <fstream>
#include "../surchem.h"
#include "../../module_xc/xc_functional.h"
#include "setcell.h"
/************************************************
 *  unit test of functions in cal_vcav.cpp
 ***********************************************/

/**
 * - Tested functions in cal_vcav.cpp:
 *   - lapl_rho
 *     - calculate Laplacian(rho)
 *   - shape_gradn
 *     - calculate first derivative of the shape function
 *   - createcavity
 *     - calculate cavitation energy
 *   - cal_vcav
 *     - calculate cavitation potential
 */

class cal_vcav_test : public testing::Test
{
protected:
    
    surchem solvent_model;
};
TEST_F(cal_vcav_test, lapl_rho)
{   
    Setcell::setupcell(GlobalC::ucell);

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
    GlobalC::rhopw ->initgrids(GlobalC::ucell.lat0,GlobalC::ucell.latvec,wfcecut);

    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();
    GlobalC::rhopw ->collect_uniqgg();

    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;
    std::complex<double> *gdrtmpg = new std::complex<double>[npw];
    double *lapn = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(lapn, nrxx);

    std::complex<double> *aux = new std::complex<double>[GlobalC::rhopw->nmaxgr];

    gdrtmpg[0] = {2.431e-07,4.760e-08};
    gdrtmpg[1] = {-7.335e-08,9.826e-07};
    gdrtmpg[2] = {-5.418e-06,2.577e-07};

    for (int ig = 3; ig < npw; ig++)
    {
        gdrtmpg[ig] = 0;
    }
    for(int i = 0 ; i < 3 ; ++i)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            aux[ig] = gdrtmpg[ig] * pow(GlobalC::rhopw->gcar[ig][i], 2);
        }
        GlobalC::rhopw->recip2real(aux, aux);
        for (int ir = 0; ir < nrxx; ir++)
        {
            lapn[ir] -= aux[ir].real() * GlobalC::ucell.tpiba2;
        }
    }

    EXPECT_NEAR(lapn[0],0.0002940907,1e-10);
    EXPECT_NEAR(lapn[1],-0.000271296,1e-10);

    delete[] gdrtmpg;
    delete[] aux;
}

TEST_F(cal_vcav_test, shape_gradn)
{  
    int nrxx = 27000;
    double TWO_PI = 6.283;
    double sigma_k = 0.6;
    double nc_k = 3.7e-04;
    double *PS_TOTN_real = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN_real,nrxx);

    PS_TOTN_real[0] = 2.081e-03; PS_TOTN_real[1] = 1.818e-03; PS_TOTN_real[2] = 1.193e-03;
    for (int i = 3; i < nrxx; i++)
    {
        PS_TOTN_real[i] = 0.1; 
    }

    double epr_c = 1.0 / sqrt(TWO_PI) / sigma_k;
    double epr_z = 0;
    double min = 1e-10;
    double *eprime = new double[nrxx];

    for (int ir = 0; ir < nrxx; ir++)
    {
        epr_z = log(max(PS_TOTN_real[ir], min) / nc_k) / sqrt(2) / sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / max(PS_TOTN_real[ir], min);
    }

    EXPECT_NEAR(eprime[0],5.0729550913,1e-10);
    EXPECT_NEAR(eprime[1],10.8252209597,1e-10);

    delete[] PS_TOTN_real;
}

TEST_F(cal_vcav_test, createcavity)
{   
    Setcell::setupcell(GlobalC::ucell);

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
    GlobalC::rhopw ->initgrids(GlobalC::ucell.lat0,GlobalC::ucell.latvec,wfcecut);

    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();
    GlobalC::rhopw ->collect_uniqgg();

    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;
    std::complex<double> *PS_TOTN = new std::complex<double>[npw];
    double *vwork = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(vwork, nrxx);

    PS_TOTN[0] = {2.432e-07,4.862e-08};
    PS_TOTN[1] = {-7.649e-08,9.806e-07};
    PS_TOTN[2] = {-5.407e-06,2.549e-07};

    for (int ig = 3; ig < npw; ig++)
    {
        PS_TOTN[ig] = 1e-7;
    }
    
    solvent_model.createcavity(GlobalC::ucell,GlobalC::rhopw,PS_TOTN,vwork);

    EXPECT_NEAR(vwork[0], 4.8556305312,1e-10);
    EXPECT_NEAR(vwork[1],-2.1006480538,1e-10);

    delete[] PS_TOTN;
    delete[] vwork;
}

TEST_F(cal_vcav_test, cal_vcav)
{   
    Setcell::setupcell(GlobalC::ucell);

    std::string precision_flag, device_flag;
    precision_flag = "double";
    device_flag = "cpu";

    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    GlobalC::rhopw = &pwtest;
    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  
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
    GlobalC::rhopw ->initgrids(GlobalC::ucell.lat0,GlobalC::ucell.latvec,wfcecut);

    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();
    GlobalC::rhopw ->collect_uniqgg();

    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;
    std::complex<double> *PS_TOTN = new std::complex<double>[npw];

    PS_TOTN[0] = {2.432e-07,4.862e-08};
    PS_TOTN[1] = {-7.649e-08,9.806e-07};
    PS_TOTN[2] = {-5.407e-06,2.549e-07};

    for (int ig = 3; ig < npw; ig++)
    {
        PS_TOTN[ig] = 1e-7;
    }
    
    int nspin = 2;
    solvent_model.Vcav.create(nspin, nrxx);

    solvent_model.cal_vcav(GlobalC::ucell,GlobalC::rhopw,PS_TOTN,nspin);

    EXPECT_NEAR(solvent_model.Vcav(0,0), 4.8556305312,1e-10);
    EXPECT_NEAR(solvent_model.Vcav(0,1),-2.1006480538,1e-10);
    EXPECT_NEAR(solvent_model.Vcav(1,0), 4.8556305312,1e-10);
    EXPECT_NEAR(solvent_model.Vcav(1,1),-2.1006480538,1e-10);

    delete[] PS_TOTN;
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
