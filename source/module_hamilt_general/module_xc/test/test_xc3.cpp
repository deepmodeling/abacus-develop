#include "../xc_functional.h"
#include "gtest/gtest.h"
#include "../exx_info.h"
#include "xc3_mock.h"
#include "module_base/matrix.h"

class XCTest_GRADCORR : public testing::Test
{
    protected:
        std::vector<double> e_lda, v_lda;
        std::vector<double> e_gga, v1_gga, v2_gga;

        void SetUp()
        {
            ModulePW::PW_Basis rhopw;
            UnitCell ucell;
            Charge chr;

            double et, vt;
            ModuleBase::matrix v;
            std::vector<double> stress;

            XC_Functional::gradcorr(et,vt,v,&chr,&rhopw,&ucell,stress,false);
        }
};

TEST_F(XCTest_GRADCORR, set_xc_type)
{

}
