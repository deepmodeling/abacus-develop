#include "../formatter_contextfmt.h"
#include <gtest/gtest.h>

/************************************************
 *  unit test of class ContextFmt
 ***********************************************/

/**
 * - Tested Functions:
 *   - ContextFmt()
 *     - Default constructor
 *   - set_context()
 *     - Set context, adjust PhysicalFmt in it
 *   - operator<<()
 *     - Stream operator, import data, convert format and save data into Table object
 *   - reset()
 *     - Reset the context
 *   - str()
 *     - Pop up the data in Table object imported by using << operator
 *   - enable_title()
 *     - Enable title of Table object
 *   - disable_title()
 *     - Disable title of Table object
 *   - only_title()
 *     - Only print title of Table object
 *   - enable_iterative()
 *     - For SCF case only presently, will print title for once and data for many times
 *   - disable_iterative()
 *     - Disable iterative
 */

class ContextFmtTest : public testing::Test
{
};

TEST_F(ContextFmtTest, DefaultConstructor) {
    formatter::ContextFmt context_fmt;
    EXPECT_NE(context_fmt.get_p_phys_fmt(), nullptr);
    EXPECT_EQ(context_fmt.get_p_phys_fmt()->get_context(), "none");
    EXPECT_EQ(context_fmt.get_mode(), 1);
    EXPECT_EQ(context_fmt.get_ncol(), 0);
    EXPECT_EQ(context_fmt.get_phys_fmt().size(), 0);
    EXPECT_EQ(context_fmt.get_nrows().size(), 0);
}

TEST_F(ContextFmtTest, SetContext1branch1) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("test");
    EXPECT_EQ(context_fmt.get_p_phys_fmt()->get_context(), "none");
    EXPECT_EQ(context_fmt.get_ncol(), 0);
    EXPECT_EQ(context_fmt.get_phys_fmt().size(), 0);
}

TEST_F(ContextFmtTest, SectContext1branch2) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    EXPECT_EQ(context_fmt.get_context(), "vector3d");
    EXPECT_EQ(context_fmt.get_ncol(), 3);
    EXPECT_EQ(context_fmt.get_phys_fmt().size(), 3);
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(context_fmt.get_phys_fmt()[i], "coordinate");
    }
    EXPECT_EQ(context_fmt.get_mode(), 1);
}

TEST_F(ContextFmtTest, SetContext2) {
    formatter::ContextFmt context_fmt;
    std::string* physfmt = new std::string[3];
    physfmt[0] = "coordinate";
    physfmt[1] = "coordinate";
    physfmt[2] = "coordinate";
    int* nrow = new int[3];
    nrow[0] = 3;
    nrow[1] = 2;
    nrow[2] = 4;
    context_fmt.set_context(3, physfmt, nrow);
    EXPECT_EQ(context_fmt.get_context(), "none");
    EXPECT_EQ(context_fmt.get_ncol(), 3);
    EXPECT_EQ(context_fmt.get_phys_fmt().size(), 3);
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(context_fmt.get_phys_fmt()[i], "coordinate");
    }
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(context_fmt.get_nrows()[i], nrow[i]);
    }
    delete[] physfmt;
    delete[] nrow;
}

TEST_F(ContextFmtTest, SetContext3) {
    formatter::ContextFmt context_fmt;
    std::vector<std::string> physfmt;
    physfmt.push_back("coordinate");
    physfmt.push_back("coordinate");
    physfmt.push_back("coordinate");
    context_fmt.set_context(physfmt);
    EXPECT_EQ(context_fmt.get_context(), "none");
    EXPECT_EQ(context_fmt.get_ncol(), 3);
    EXPECT_EQ(context_fmt.get_phys_fmt().size(), 3);
    for (int i = 0; i < 3; ++i) {
        EXPECT_EQ(context_fmt.get_phys_fmt()[i], "coordinate");
    }
}
// vector3d context
TEST_F(ContextFmtTest, StreamOperatorLeft1int) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    context_fmt<<2<<3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft1double) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<1.0;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    context_fmt<<2.0<<3.0;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft1string) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<"1";
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    context_fmt<<"2"<<"3";
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft1char) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<'1';
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    context_fmt<<'2'<<'3';
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2VecInt) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    std::vector<int> v1 = {1, 2, 3};
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::vector<int> v2 = {4, 5, 6};
    std::vector<int> v3 = {7, 8, 9};
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2VecDouble) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    std::vector<double> v1 = {1.0, 2.0, 3.0};
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::vector<double> v2 = {4.0, 5.0, 6.0};
    std::vector<double> v3 = {7.0, 8.0, 9.0};
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2VecString) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    std::vector<std::string> v1 = {"1", "2", "3"};
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::vector<std::string> v2 = {"4", "5", "6"};
    std::vector<std::string> v3 = {"7", "8", "9"};
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2PtrInt) {
    formatter::ContextFmt context_fmt;
    int* v_nrows = new int[3];
    v_nrows[0] = 3;
    v_nrows[1] = 3;
    v_nrows[2] = 3;
    context_fmt.set_context("vector3d", 3, v_nrows);
    int* v1 = new int[3];
    v1[0] = 1;
    v1[1] = 2;
    v1[2] = 3;
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    int* v2 = new int[3];
    v2[0] = 4;
    v2[1] = 5;
    v2[2] = 6;
    int* v3 = new int[3];
    v3[0] = 7;
    v3[1] = 8;
    v3[2] = 9;
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    delete[] v_nrows;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2PtrDouble) {
    formatter::ContextFmt context_fmt;
    int* v_nrows = new int[3];
    v_nrows[0] = 3;
    v_nrows[1] = 3;
    v_nrows[2] = 3;
    context_fmt.set_context("vector3d", 3, v_nrows);
    double* v1 = new double[3];
    v1[0] = 1.;
    v1[1] = 2.;
    v1[2] = 3.;
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    double* v2 = new double[3];
    v2[0] = 4.;
    v2[1] = 5.;
    v2[2] = 6.;
    double* v3 = new double[3];
    v3[0] = 7.;
    v3[1] = 8.;
    v3[2] = 9.;
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    delete[] v_nrows;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2PtrString) {
    formatter::ContextFmt context_fmt;
    int* v_nrows = new int[3];
    v_nrows[0] = 3;
    v_nrows[1] = 3;
    v_nrows[2] = 3;
    context_fmt.set_context("vector3d", 3, v_nrows);
    std::string* v1 = new std::string[3];
    v1[0] = "1";
    v1[1] = "2";
    v1[2] = "3";
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    std::string* v2 = new std::string[3];
    v2[0] = "4";
    v2[1] = "5";
    v2[2] = "6";
    std::string* v3 = new std::string[3];
    v3[0] = "7";
    v3[1] = "8";
    v3[2] = "9";
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    delete[] v_nrows;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    std::cout<<context_fmt.str()<<std::endl;
}

TEST_F(ContextFmtTest, StreamOperatorLeft2PtrChar) {
    formatter::ContextFmt context_fmt;
    int* v_nrows = new int[3];
    v_nrows[0] = 3;
    v_nrows[1] = 3;
    v_nrows[2] = 3;
    context_fmt.set_context("vector3d", 3, v_nrows);
    char* v1 = new char[3];
    v1[0] = '1';
    v1[1] = '2';
    v1[2] = '3';
    context_fmt<<v1;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 1);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    char* v2 = new char[3];
    v2[0] = '4';
    v2[1] = '5';
    v2[2] = '6';
    char* v3 = new char[3];
    v3[0] = '7';
    v3[1] = '8';
    v3[2] = '9';
    context_fmt<<v2<<v3;
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 3);
    EXPECT_EQ(context_fmt.get_title_switch()%2, 1);
    delete[] v_nrows;
    delete[] v1;
    delete[] v2;
    delete[] v3;
    std::cout<<context_fmt.str()<<std::endl;
    std::cout<<"Note: ASCII code of '1' is "<<(int)'1'<<std::endl;
}

TEST_F(ContextFmtTest, Reset) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<1<<2<<3;
    context_fmt.reset();
    EXPECT_EQ(context_fmt.get_cache_title(), "");
    EXPECT_EQ(context_fmt.get_icol(), 0);
    EXPECT_EQ(context_fmt.get_title_switch(), 1);
}

TEST_F(ContextFmtTest, Str) {
    formatter::ContextFmt context_fmt;
    context_fmt.set_context("vector3d");
    context_fmt<<1<<2<<3;
    EXPECT_EQ(context_fmt.str(), "    1.0000000000     2.0000000000     3.0000000000");
}

TEST_F(ContextFmtTest, EnableTitle) {
    formatter::ContextFmt context_fmt;
    context_fmt.enable_title();
    EXPECT_EQ(context_fmt.get_mode(), 0);
    EXPECT_EQ(context_fmt.get_title_switch(), 0);
}

TEST_F(ContextFmtTest, DisableTitle) {
    formatter::ContextFmt context_fmt;
    context_fmt.disable_title();
    EXPECT_EQ(context_fmt.get_mode(), 1);
    EXPECT_EQ(context_fmt.get_title_switch(), 1);
}

TEST_F(ContextFmtTest, OnlyTitle) {
    formatter::ContextFmt context_fmt;
    context_fmt.only_title();
    EXPECT_EQ(context_fmt.get_mode(), -1);
    EXPECT_EQ(context_fmt.get_title_switch(), 0);
}

TEST_F(ContextFmtTest, EnableIterative) {
    formatter::ContextFmt context_fmt;
    context_fmt.enable_iterative();
    EXPECT_EQ(context_fmt.get_iterative(), 1);
}

TEST_F(ContextFmtTest, DisableIterative) {
    formatter::ContextFmt context_fmt;
    context_fmt.disable_iterative();
    EXPECT_EQ(context_fmt.get_iterative(), 0);
}
