#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>
#include <fstream>
#include "../surchem.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"

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
    
    double *PS_TOTN_real = new double[GlobalC::rhopw->nrxx];
    surchem solvent_model;
    
};

TEST_F(cal_epsilon_test, cal_epsilon)
{   
    ifstream fin;
	fin.open("./support/PS_TOTN_real");
    for (int i=0; i<GlobalC::rhopw->nrxx; i++)
    {
        fin>>PS_TOTN_real[i];
    }
    GlobalC::rhopw->nrxx = 5;
    double *epsilon = new double[GlobalC::rhopw->nrxx];
    double *epsilon0 = new double[GlobalC::rhopw->nrxx];
    solvent_model.cal_epsilon(GlobalC::rhopw, PS_TOTN_real, epsilon, epsilon0);
    EXPECT_EQ(epsilon[0], 1);
    EXPECT_EQ(epsilon[12],1.00005);
    EXPECT_EQ(epsilon[19], 43.1009);
    EXPECT_EQ(epsilon[26], 78.746);
}