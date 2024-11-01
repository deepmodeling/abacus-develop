#include <gtest/gtest.h>
#include "module_hamilt_general/module_vdw/dftd3_xc_name.h"
#include "module_hamilt_general/module_vdw/dftd3_xc_param.h"

TEST(DFTD3XCTest, SearchXcnameLibXCXplusC)
{
    std::string xcpattern = "XC_GGA_X_PBE+XC_GGA_C_OP_PBE";
    std::string xname;
    DFTD3::search_xcname_libxc_xplusc(xcpattern, xname);
    EXPECT_EQ(xname, "pbeop");
}

TEST(DFTD3XCTest, SearchXcnameLibXCXC)
{
    std::string xcpattern = "XC_LDA_XC_TETER93";
    std::string xname;
    DFTD3::search_xcname_libxc_xc(xcpattern, xname);
    EXPECT_EQ(xname, "teter93");
}

TEST(DFTD3XCTest, SearchXcnameLibXC)
{
    std::string xcpattern = "XC_GGA_X_PBE+XC_GGA_C_OP_PBE";
    std::string xname;
    DFTD3::search_xcname_libxc(xcpattern, xname);
    EXPECT_EQ(xname, "pbeop");
}

TEST(DFTD3XCTest, SearchXcname)
{
    std::string xcpattern = "XC_GGA_X_PBE+XC_GGA_C_OP_PBE";
    std::string xname = DFTD3::search_xcname(xcpattern);
    EXPECT_EQ(xname, "pbeop");

    xcpattern = "XC_LDA_XC_TETER93";
    xname = DFTD3::search_xcname(xcpattern);
    EXPECT_EQ(xname, "teter93");

    xcpattern = "default";
    xname = DFTD3::search_xcname(xcpattern);
    EXPECT_EQ(xname, "default");

    xcpattern = "pbe";
    xname = DFTD3::search_xcname(xcpattern);
    EXPECT_EQ(xname, "pbe");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}