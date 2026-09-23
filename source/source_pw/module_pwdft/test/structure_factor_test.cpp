#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <string>
#include <cmath>
#include <complex>
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/test/prepare_unitcell.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/stru_fac.h"
/************************************************
 *  unit test of class Structure_factor and 
 ***********************************************/

/**
 * - Tested Functions:
 *   - Fcoef::create to create a 5 dimensional array of complex numbers
 *   - Soc::set_fcoef to set the fcoef array
 *   - Soc::spinor to calculate the spinor
 *   - Soc::rot_ylm to calculate the rotation matrix
 *   - Soc::sph_ind to calculate the m index of the spherical harmonics
*/

//compare two complex by using EXPECT_DOUBLE_EQ()


Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

class StructureFactorTest : public testing::Test
{
protected:
    Structure_Factor SF;
    std::string output;
    ModulePW::PW_Basis* rho_basis;
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    Parallel_Grid* pgrid;
    std::vector<int> nw = {13};
    int nlocal = 0;
void SetUp()
{
    rho_basis=new ModulePW::PW_Basis;
    ucell = utp.SetUcellInfo(nw, nlocal);
    ucell->set_iat2iwt(1);
    pgrid = new Parallel_Grid;
    rho_basis->npw=10;
    rho_basis->gcar=new ModuleBase::Vector3<double>[10];
    // for (int ig=0;ig<rho_basis->npw;ig++)
    // {
    //     rho_basis->gcar[ig]=1.0;
    // }
}
};

TEST_F(StructureFactorTest, set)
{
    const ModulePW::PW_Basis* rho_basis_in;
    const int nbspline_in =10;
    SF.set(rho_basis_in,nbspline_in);
    EXPECT_EQ(nbspline_in, 10);
}


TEST_F(StructureFactorTest, setup_structure_factor_double)
{
    rho_basis->npw = 10;
    SF.setup(ucell,*pgrid,rho_basis,false);

    for (int i=0;i< ucell->nat * (2 * rho_basis->nx + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts1_data<double>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts1_data<double>()[i].imag(),0);
    }

    for (int i=0;i< ucell->nat * (2 * rho_basis->ny + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts2_data<double>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts2_data<double>()[i].imag(),0);
    }

    for (int i=0;i< ucell->nat * (2 * rho_basis->nz + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts3_data<double>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts3_data<double>()[i].imag(),0);
    }
}

TEST_F(StructureFactorTest, setup_structure_factor_float)
{
    // the float eigts copies are what this case checks, so ask setup() for them
    rho_basis->npw = 10;
    SF.setup(ucell,*pgrid,rho_basis,true);

    for (int i=0;i< ucell->nat * (2 * rho_basis->nx + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts1_data<float>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts1_data<float>()[i].imag(),0);
    }

    for (int i=0;i< ucell->nat * (2 * rho_basis->ny + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts2_data<float>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts2_data<float>()[i].imag(),0);
    }

    for (int i=0;i< ucell->nat * (2 * rho_basis->nz + 1);i++) 
    {
       EXPECT_EQ(SF.get_eigts3_data<float>()[i].real(),1);
       EXPECT_EQ(SF.get_eigts3_data<float>()[i].imag(),0);
    }
}

int main()
{
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
