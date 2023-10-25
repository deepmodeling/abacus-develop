#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

K_Vectors::K_Vectors(){}
K_Vectors::~K_Vectors(){}

/************************************************
 *  unit test of the function of cal_h_lambda
 ***********************************************/

/**
 * Tested function:
 *  - SpinConstrain::cal_h_lambda
 *    - this function calculates the h_lambda operator
 */

class SpinConstrainTest : public testing::Test
{
  protected:
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, CalHLambda)
{
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    int nrow = 2;
    int ncol = 2;
    paraV.set_global2local(nrow, ncol, false, ofs);
    sc.set_ParaV(&paraV);
    EXPECT_EQ(sc.ParaV->nloc, 4);
    std::map<int, int> atomCounts = {
        {0, 1}
    };
    std::map<int, int> orbitalCounts = {
        {0, 1}
    };
    sc.clear_atomCounts();
    sc.clear_orbitalCounts();
    sc.set_atomCounts(atomCounts);
    sc.set_orbitalCounts(orbitalCounts);
    sc.set_npol(2);
    std::vector<std::complex<double>> h_lambda(sc.ParaV->nloc);
    std::fill(h_lambda.begin(), h_lambda.end(), std::complex<double>(0, 0));
    std::vector<std::complex<double>> Sloc2 = {
        std::complex<double>{1.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{0.0, 0.0},
        std::complex<double>{1.0, 0.0}
    };
    ModuleBase::Vector3<double>* sc_lambda = new ModuleBase::Vector3<double>[1];
    sc_lambda[0][0] = 1.0;
    sc_lambda[0][1] = 1.0;
    sc_lambda[0][2] = 1.0;
    sc.set_sc_lambda(sc_lambda, 1);
    sc.cal_h_lambda(&h_lambda[0], Sloc2, true);
    delete[] sc_lambda;
    remove("test.log");
}