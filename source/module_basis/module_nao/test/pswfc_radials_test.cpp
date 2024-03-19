#include <gtest/gtest.h>
#include "module_basis/module_nao/pswfc_radials.h"
#include <fstream>

#ifdef __MPI
#include <mpi.h>
#endif

class PswfcRadialsTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
#ifdef __MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif
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

    pswfc_radials.read_upf_pswfc(ifs, 0.0, 1e-6);
    EXPECT_EQ(pswfc_radials.lmax(), 2);
    EXPECT_EQ(pswfc_radials.nzeta(0), 1);
    EXPECT_EQ(pswfc_radials.nzeta(1), 1);
    EXPECT_EQ(pswfc_radials.nzeta(2), 1);
    EXPECT_EQ(pswfc_radials.nzeta_max(), 1);
    // l = 0, <PP_CHI.2>, 4S
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(0), 5.0672226831E-13);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(1), 3.0740550920E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(2), 6.2055866358E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(3), 9.4519832136E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(4), 1.2870457911E-03);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(39), 6.2747938942E-02);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(79), 3.2957297188E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(119), 6.5729325723E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(159), 7.9744230720E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(0, 0).rvalue(199), 7.3512466100E-01);
    // l = 1, <PP_CHI.3>, 4P
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(0), -1.2105620326E-12);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(1), 2.6035363423E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(2), 1.0416837423E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(3), 2.3447967866E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(4), 4.1710331234E-04);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(39), 4.3850721200E-02);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(79), 2.0385391536E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(119), 4.4174606396E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(159), 6.1861552438E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(1, 0).rvalue(199), 6.6885272981E-01);
    // l = 2, <PP_CHI.1>, 3D
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(0), 2.0998757420E-11);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(1), 1.0346699413E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(2), 8.2719900770E-05);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(3), 2.7887785682E-04);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(4), 6.6004347176E-04);

    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(39), 4.4075591655E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(79), 1.2790898686E+00);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(119), 7.5247184647E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(159), 2.7836191101E-01);
    EXPECT_DOUBLE_EQ(pswfc_radials.chi(2, 0).rvalue(199), 1.3889295980E-01);
}

TEST_F(PswfcRadialsTest, Build)
{
    PswfcRadials pswfc_radials;
    std::string pspot = "../../../../../tests/PP_ORB/As_dojo.upf";
    pswfc_radials.build(pspot, 0, 0.0);
    
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

#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    std::cout << "Current getcwd: " << getcwd(nullptr, 0) << std::endl;
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}