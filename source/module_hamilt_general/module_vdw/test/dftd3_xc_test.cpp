#include <gtest/gtest.h>
#include "module_hamilt_general/module_vdw/dftd3_xc_name.h"
#include "module_hamilt_general/module_vdw/dftd3_xc_param.h"

TEST(DFTD3XCTest, SearchXcnameLibXCXplusC)
{
    std::string xname;
    DFTD3::_xcname_libxc_xplusc("XC_GGA_X_PBE+XC_GGA_C_OP_PBE", xname);
    EXPECT_EQ(xname, "pbeop");
    // then test the case with out prefix XC_
    DFTD3::_xcname_libxc_xplusc("GGA_X_PBE+GGA_C_OP_PBE", xname);
    EXPECT_EQ(xname, "pbeop");
}

TEST(DFTD3XCTest, SearchXcnameLibXCXC)
{
    std::string xname;
    DFTD3::_xcname_libxc_xc("XC_LDA_XC_TETER93", xname);
    EXPECT_EQ(xname, "teter93");
    // then test the case with out prefix XC_
    DFTD3::_xcname_libxc_xc("LDA_XC_TETER93", xname);
    EXPECT_EQ(xname, "teter93");
}

TEST(DFTD3XCTest, SearchXcnameLibXC)
{
    std::string xname;
    DFTD3::_xcname_libxc("XC_GGA_X_PBE+XC_GGA_C_OP_PBE", xname);
    EXPECT_EQ(xname, "pbeop");
    // then test the case with out prefix XC_
    DFTD3::_xcname_libxc("GGA_X_PBE+GGA_C_OP_PBE", xname);
    EXPECT_EQ(xname, "pbeop");
}

TEST(DFTD3XCTest, SearchXcname)
{
    std::string xcpattern = "XC_GGA_X_PBE+XC_GGA_C_OP_PBE";
    std::string xname = DFTD3::_xcname(xcpattern);
    EXPECT_EQ(xname, "pbeop");

    xcpattern = "XC_LDA_XC_TETER93";
    xname = DFTD3::_xcname(xcpattern);
    EXPECT_EQ(xname, "teter93");

    xcpattern = "default";
    xname = DFTD3::_xcname(xcpattern);
    EXPECT_EQ(xname, "default");

    xcpattern = "pbe";
    xname = DFTD3::_xcname(xcpattern);
    EXPECT_EQ(xname, "pbe");
}

TEST(DFTD3XCTest, SuccessfulSearch)
{
    std::string xc = "pbe";
    std::string d3method = "d3_0";
    std::string s6_in = "default";
    std::string s8_in = "default";
    std::string a1_in = "default";
    std::string a2_in = "default";
    double s6, s8, a1, a2;
    DFTD3::dftd3_params(xc, d3method, s6_in, s8_in, a1_in, a2_in, s6, s8, a1, a2);
    EXPECT_DOUBLE_EQ(s6, 1.0);
    EXPECT_DOUBLE_EQ(s8, 0.722);
    EXPECT_DOUBLE_EQ(a1, 1.217);
    EXPECT_DOUBLE_EQ(a2, 1.0);

    // a more complicated case: MGGA_X_SCAN+MGGA_C_SCAN
    xc = "XC_MGGA_X_SCAN+XC_MGGA_C_SCAN";
    DFTD3::dftd3_params(xc, d3method, s6_in, s8_in, a1_in, a2_in, s6, s8, a1, a2);
    EXPECT_DOUBLE_EQ(s6, 1.0);
    EXPECT_DOUBLE_EQ(s8, 0.0);
    EXPECT_DOUBLE_EQ(a1, 1.324);
    EXPECT_DOUBLE_EQ(a2, 1.0);

    // user defines all parameters
    s6_in = "1.1";
    s8_in = "0.1";
    a1_in = "1.325";
    a2_in = "1.1";
    DFTD3::dftd3_params(xc, d3method, s6_in, s8_in, a1_in, a2_in, s6, s8, a1, a2);
    EXPECT_DOUBLE_EQ(s6, 1.1);
    EXPECT_DOUBLE_EQ(s8, 0.1);
    EXPECT_DOUBLE_EQ(a1, 1.325);
    EXPECT_DOUBLE_EQ(a2, 1.1);

    // user defines one parameter
    s6_in = "1.1";
    s8_in = "default";
    a1_in = "default";
    a2_in = "default";
    DFTD3::dftd3_params(xc, d3method, s6_in, s8_in, a1_in, a2_in, s6, s8, a1, a2);
    EXPECT_DOUBLE_EQ(s6, 1.1);
    EXPECT_DOUBLE_EQ(s8, 0.0);
    EXPECT_DOUBLE_EQ(a1, 1.324);
    EXPECT_DOUBLE_EQ(a2, 1.0);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}