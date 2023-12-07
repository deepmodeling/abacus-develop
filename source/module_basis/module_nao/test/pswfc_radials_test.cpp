#include <gtest/gtest.h>
#include "module_basis/module_nao/pswfc_radials.h"
#include <fstream>

class PswfcRadialsTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
            // code here will execute just before the test ensues 
            // (each TEST_F)
        }
        virtual void TearDown() {
            // code here will be called just after the test completes
            // (each TEST_F)
        }
};

TEST_F(PswfcRadialsTest, startswith) {
    PswfcRadials pswfc_radials;
    std::string word = "hello";
    std::string pattern = "he";
    bool result = pswfc_radials.startswith(word, pattern);
    EXPECT_TRUE(result);
}

TEST_F(PswfcRadialsTest, StealFromQuotes) {
    PswfcRadials pswfc_radials;
    std::string word = "\"hello\"";
    std::string result = pswfc_radials.steal_from_quotes(word);
    EXPECT_EQ(result, "hello");
}

TEST_F(PswfcRadialsTest, StealFromQuotesOverload1) {
    PswfcRadials pswfc_radials;
    std::string pspot = "../../../../../tests/PP_ORB/As_dojo.upf";
    std::ifstream ifs;
    ifs.open(pspot);
    // check if file is opened
    bool is_open = ifs.is_open();
    if(!is_open) {std::cout<<"File path WRONG.\n"; }
    ASSERT_TRUE(is_open);

    std::string line;
    while(!ifs.eof())
    {
        ifs >> line;
        if(pswfc_radials.startswith(line, "mesh_size="))
        {
            int result = std::stoi(pswfc_radials.steal_from_quotes(ifs, line));
            EXPECT_EQ(result, 1358);
        }
    }
    ifs.close();
}

TEST_F(PswfcRadialsTest, ReadKeywordValue) {
    PswfcRadials pswfc_radials;
    std::string pspot = "../../../../../tests/PP_ORB/As_dojo.upf";
    std::ifstream ifs;
    ifs.open(pspot);
    // check if file is opened
    bool is_open = ifs.is_open();
    if(!is_open) {std::cout<<"File path WRONG.\n"; }
    ASSERT_TRUE(is_open);

    std::string line;
    while(!ifs.eof())
    {
        ifs >> line;
        if(pswfc_radials.startswith(line, "mesh_size="))
        {
            std::string str_result = pswfc_radials.read_keyword_value(ifs, line);
            int result = std::stoi(str_result);
            EXPECT_EQ(result, 1358);
        }
    }
    ifs.close();
}

TEST_F(PswfcRadialsTest, ReadUpfPswfc) {
    PswfcRadials pswfc_radials;
    std::string pspot = "../../../../../tests/PP_ORB/As_dojo.upf";
    std::ifstream ifs;
    ifs.open(pspot);
    // check if file is opened
    bool is_open = ifs.is_open();
    if(!is_open) {std::cout<<"File path WRONG.\n"; }
    ASSERT_TRUE(is_open);

    pswfc_radials.read_upf_pswfc(ifs);
    EXPECT_EQ(pswfc_radials.lmax(), 2);
    EXPECT_EQ(pswfc_radials.nzeta(0), 1);
    EXPECT_EQ(pswfc_radials.nzeta(1), 1);
    EXPECT_EQ(pswfc_radials.nzeta(2), 1);
    EXPECT_EQ(pswfc_radials.nzeta_max(), 1);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(0), 5.0672226831E-13);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(1), 3.0740550920E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(2), 6.2055866358E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(3), 9.4519832136E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(4), 1.2870457911E-03);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(0), -1.2105620326E-12);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(1), 2.6035363423E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(2), 1.0416837423E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(3), 2.3447967866E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(4), 4.1710331234E-04);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(0), 2.0998757420E-11);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(1), 1.0346699413E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(2), 8.2719900770E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(3), 2.7887785682E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(4), 6.6004347176E-04);
}

TEST_F(PswfcRadialsTest, Build)
{
    PswfcRadials pswfc_radials;
    std::string pspot = "../../../../../tests/PP_ORB/As_dojo.upf";
    pswfc_radials.build(pspot);
    
    EXPECT_EQ(pswfc_radials.lmax(), 2);
    EXPECT_EQ(pswfc_radials.nzeta(0), 1);
    EXPECT_EQ(pswfc_radials.nzeta(1), 1);
    EXPECT_EQ(pswfc_radials.nzeta(2), 1);
    EXPECT_EQ(pswfc_radials.nzeta_max(), 1);
    EXPECT_EQ(pswfc_radials.rcut_max(), 13.57);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(0), 5.0672226831E-13);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(1), 3.0740550920E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(2), 6.2055866358E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(3), 9.4519832136E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(4), 1.2870457911E-03);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(0), -1.2105620326E-12);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(1), 2.6035363423E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(2), 1.0416837423E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(3), 2.3447967866E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(4), 4.1710331234E-04);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(0), 2.0998757420E-11);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(1), 1.0346699413E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(2), 8.2719900770E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(3), 2.7887785682E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(4), 6.6004347176E-04);
}

int main(int argc, char** argv)
{
    std::cout << "Current getcwd: " << getcwd(nullptr, 0) << std::endl;
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    return result;
}