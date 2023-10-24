#include "../spin_constrain.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "prepare_unitcell.h"

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

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
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    // nw is the number of orbitals of each atom
    // it should container ucell.nat elements
    std::vector<int> nw = {13};
    int nlocal;
    void SetUp() override
    {
        ucell = utp.SetUcellInfo(nw, nlocal);
    }
    // get sc instance
    SpinConstrain<std::complex<double>, psi::DEVICE_CPU>& sc
        = SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::getScInstance();
};

TEST_F(SpinConstrainTest, Foo)
{
}