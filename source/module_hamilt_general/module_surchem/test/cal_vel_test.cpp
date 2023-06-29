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
 *  unit test of functions in cal_vel.cpp
 ***********************************************/

/**
 * - Tested functions in cal_vel.cpp:
 *   - shape_gradn
 *     - calculate first derivative of the shape function
 *   - eps_pot
 *     - calculate the 2nd item of Vel
 *   - cal_vel
 *     - calculate electrostatic potential
 */

class cal_vel_test : public testing::Test
{
protected:
    
    surchem solvent_model;
};

TEST_F(cal_vel_test, shape_gradn)
{  
    int nrxx = 27000;
    double TWO_PI = 6.283;
    double sigma_k = 0.6;
    double nc_k = 3.7e-04;
    double *PS_TOTN_real = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN_real,nrxx);

    PS_TOTN_real[0]=2.081e-03; PS_TOTN_real[1]=1.818e-03;PS_TOTN_real[2]=1.193e-03;
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
        epr_z = log(max(PS_TOTN_real[ir], min) / GlobalV::nc_k) / sqrt(2) / GlobalV::sigma_k;
        eprime[ir] = epr_c * exp(-pow(epr_z, 2)) / max(PS_TOTN_real[ir], min);
    }

    EXPECT_NEAR(eprime[0],5.0729550913,1e-10);
    EXPECT_NEAR(eprime[1],10.8252209597,1e-10);

    delete[] PS_TOTN_real;
}

TEST_F(cal_vel_test, eps_pot)
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
    double *vwork = new double[nrxx];
    double *eprime = new double[nrxx];
    double *PS_TOTN_real = new double[nrxx];
    ModuleBase::GlobalFunc::ZEROS(vwork, nrxx);
    ModuleBase::GlobalFunc::ZEROS(eprime, nrxx);
    ModuleBase::GlobalFunc::ZEROS(PS_TOTN_real,nrxx);
    
    PS_TOTN_real[0] = 2.081e-03; PS_TOTN_real[1] = 1.818e-03; PS_TOTN_real[2] = 1.193e-03;
    eprime[0] = 5.07288; eprime[1] = 10.8251; eprime[2] = 83.0605;
    for (int i = 3; i < nrxx; i++)
    {
        PS_TOTN_real[i] = 0.1; 
        eprime[i] = 1;
    }

    complex<double> *phi = new complex<double>[npw];
    phi[0] = {2.116e-05,-9.528e-05}; phi[1] = {1.608e-04,-1.958e-06}; phi[2] = {1.500e-05,7.303e-05};
    for(int i = 3; i < npw; i++)
    {
        phi[i]=1e-7;
    }

    for (int ir = 0; ir < nrxx; ir++)
    {
        eprime[ir] = eprime[ir] * (GlobalV::eb_k - 1);
    }

    ModuleBase::Vector3<double> *nabla_phi = new ModuleBase::Vector3<double>[nrxx];
    double *phisq = new double[nrxx];

    XC_Functional::grad_rho(phi, nabla_phi, GlobalC::rhopw, GlobalC::ucell.tpiba);

    for (int ir = 0; ir < nrxx; ir++)
    {
        phisq[ir] = pow(nabla_phi[ir].x, 2) + pow(nabla_phi[ir].y, 2) + pow(nabla_phi[ir].z, 2);
    }

    for (int ir = 0; ir < nrxx; ir++)
    {
        vwork[ir] = eprime[ir] * phisq[ir] / (8 * ModuleBase::PI);
    }

    EXPECT_NEAR(vwork[0], 1.46866e-06, 1e-10);
    EXPECT_NEAR(vwork[1], 0.0004116982, 1e-10);

    delete[] PS_TOTN_real;
    delete[] phi;
    delete[] eprime;
    delete[] nabla_phi;
    delete[] phisq;
    delete[] vwork;
}

TEST_F(cal_vel_test, cal_vel)
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

    complex<double> *TOTN = new complex<double>[npw];
    complex<double> *PS_TOTN = new complex<double>[npw];
    
    for(int i=0;i<npw;i++)
    {
        TOTN[i] = 1e-5;
        PS_TOTN[i] = 1e-7;
    }

    int nspin = 1;
    solvent_model.Vel.create(nspin, nrxx);
    solvent_model.epspot = new double[nrxx];
    solvent_model.TOTN_real =  new double[nrxx];
    solvent_model.delta_phi =  new double[nrxx];

    solvent_model.cal_vel(GlobalC::ucell, GlobalC::rhopw, TOTN, PS_TOTN, nspin);

    EXPECT_NEAR(solvent_model.Vel(0,0), 0.0532168705, 1e-10);
    EXPECT_NEAR(solvent_model.Vel(0,1), 0.0447818244, 1e-10);

    delete[] PS_TOTN;
    delete[] TOTN;
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
