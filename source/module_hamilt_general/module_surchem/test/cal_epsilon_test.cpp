#ifdef __MPI
#include "mpi.h"
#include "../../../module_basis/module_pw/test/test_tool.h"
#include "../../../module_base/parallel_global.h"
#endif
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>
#include <fstream>
#include "../surchem.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_base/constants.h"
#include "../../../module_base/parallel_global.h"

#define doublethreshold 1e-5

/************************************************
 *  unit test of functions in cal_epsilon.cpp
 ***********************************************/

/**
 * - Tested functions in cal_epsilon.cpp:
 *   - cal_epsilon
 *     - calculate the relative permittivity
 */

namespace GlobalC
{
    ModulePW::PW_Basis* rhopw;
}

class cal_epsilon_test : public testing::Test
{
protected:
    
    surchem solvent_model;
};

TEST_F(cal_epsilon_test, cal_epsilon)
{   
    std::string precision_flag, device_flag;
    precision_flag = "double";
    device_flag = "cpu";

    ModulePW::PW_Basis pwtest(device_flag, precision_flag);
    GlobalC::rhopw = &pwtest;

    ModuleBase::Matrix3 latvec;
    int nx,ny,nz;  //f*G
    double wfcecut;
    double lat0;
    bool gamma_only;
    //--------------------------------------------------
    lat0 = 1.0;
    ModuleBase::Matrix3 la(1, 1, 0, 0, 1, 1, 0, 0, 2);
    latvec = la;
    wfcecut = 60;
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
    GlobalC::rhopw ->initgrids(lat0,latvec,wfcecut);
    GlobalC::rhopw ->initparameters(gamma_only,wfcecut,distribution_type, xprime);
    GlobalC::rhopw ->setuptransform();
    GlobalC::rhopw ->collect_local_pw();

    GlobalC::rhopw ->nrxx = 125000;
    const int npw = GlobalC::rhopw ->npw;
    const int nrxx = GlobalC::rhopw ->nrxx;


    std::ifstream fin;
	fin.open("./support/PS_TOTN_real.in");
    if (!fin)
	{
	    std::cerr<<"input file does not exist"<<std::endl;
		return ;
	}
    

    double *PS_TOTN_real = new double[nrxx];
    for (int i=0; i<nrxx; i++)
    {
        fin>>PS_TOTN_real[i];
    }
    
    double *epsilon = new double[nrxx];
    double *epsilon0 = new double[nrxx];

    GlobalC::rhopw = &pwtest;
    solvent_model.cal_epsilon(GlobalC::rhopw, PS_TOTN_real, epsilon, epsilon0);

    EXPECT_EQ(PS_TOTN_real[0], 0.274231);
    EXPECT_EQ(epsilon[0],1);
    EXPECT_NEAR(epsilon[12],1.00005,doublethreshold);
    //EXPECT_EQ(epsilon[19], 43.1009);
    //EXPECT_EQ(epsilon[26], 78.746);
    delete [] PS_TOTN_real;
    delete [] epsilon;
    delete [] epsilon0;
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
