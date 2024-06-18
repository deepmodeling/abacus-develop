
#include "../output_mulliken.h"

#include "module_cell/cell_index.h"
#include "module_io/output_dmk.h"
#include "module_io/output_sk.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

template <typename TK>
class OutputMullikenTest : public testing::Test
{
  protected:
    std::vector<std::string> atomLabels = {"Si"};
    std::vector<int> atomCounts = {1};
    std::vector<std::vector<int>> lnchiCounts = {{2, 2, 1}};
    Parallel_Orbitals paraV;
    int nrow;
    int ncol;
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(OutputMullikenTest, MyTypes);

TYPED_TEST(OutputMullikenTest, OrbInfo)
{
    CellIndex cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 1);
    cell_index.write_orb_info("./");
    std::ifstream ifs("./Orbital");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("#io    spec    l    m    z  sym"));
    EXPECT_THAT(str, testing::HasSubstr("0      Si    2    3    1       dx^2-y^2"));
    remove("./Orbital");
}

#ifdef __MPI
TYPED_TEST(OutputMullikenTest, nspin1)
{
    this->nrow = 13;
    this->ncol = 13;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 1);
    auto out_sk = ModuleIO::Output_Sk<TypeParam>(nullptr, nullptr, &this->paraV, 1, 1);
    auto out_dmk = ModuleIO::Output_DMK<TypeParam>(nullptr, &this->paraV, 1, 1);
    auto mulp = ModuleIO::Output_Mulliken<TypeParam>(&(out_sk), &(out_dmk), &(this->paraV), &(cell_index), {0}, 1);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 4.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Total charge:\t4"));
    EXPECT_THAT(str, testing::HasSubstr("Total Charge on atom:                 Si              4.0000"));
    remove("./mulliken.txt");
}

TYPED_TEST(OutputMullikenTest, nspin2)
{
    this->nrow = 13;
    this->ncol = 13;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 2);
    auto out_sk = ModuleIO::Output_Sk<TypeParam>(nullptr, nullptr, &this->paraV, 2, 1);
    auto out_dmk = ModuleIO::Output_DMK<TypeParam>(nullptr, &this->paraV, 2, 1);
    auto mulp = ModuleIO::Output_Mulliken<TypeParam>(&(out_sk), &(out_dmk), &(this->paraV), &(cell_index), {0, 1}, 2);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 3.0, 1e-5);
    EXPECT_NEAR(tot_chg[1], 1.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Total charge:\t4"));
    EXPECT_THAT(str, testing::HasSubstr("Total charge of spin 1:\t3"));
    EXPECT_THAT(str, testing::HasSubstr("Total charge of spin 2:\t1"));
    EXPECT_THAT(str, testing::HasSubstr("Total Charge on atom:                 Si              4.0000"));
    EXPECT_THAT(str, testing::HasSubstr("Total Magnetism on atom:              Si              2.0000"));
    remove("./mulliken.txt");
}

TYPED_TEST(OutputMullikenTest, nspin4)
{
    this->nrow = 26;
    this->ncol = 26;
    this->paraV.init(this->nrow, this->ncol, 1, MPI_COMM_WORLD, 0);
    auto cell_index = CellIndex(this->atomLabels, this->atomCounts, this->lnchiCounts, 4);
    auto out_sk = ModuleIO::Output_Sk<std::complex<double>>(nullptr, nullptr, &this->paraV, 4, 1);
    auto out_dmk = ModuleIO::Output_DMK<std::complex<double>>(nullptr, &this->paraV, 4, 1);
    auto mulp
        = ModuleIO::Output_Mulliken<std::complex<double>>(&(out_sk), &(out_dmk), &(this->paraV), &(cell_index), {0}, 4);
    mulp.write(0, "./");
    std::vector<double> tot_chg = mulp.get_tot_chg();
    EXPECT_NEAR(tot_chg[0], 4.0, 1e-5);
    EXPECT_NEAR(tot_chg[1], 0.0, 1e-5);
    EXPECT_NEAR(tot_chg[2], 0.0, 1e-5);
    EXPECT_NEAR(tot_chg[3], 2.0, 1e-5);
    std::ifstream ifs("./mulliken.txt");
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("Total charge:\t4"));
    EXPECT_THAT(str, testing::HasSubstr("Total Charge on atom:                 Si              4.0000"));
    EXPECT_THAT(
        str,
        testing::HasSubstr(
            "Total Magnetism on atom:              Si              0.0000              0.0000              2.0000"));
    remove("./mulliken.txt");
}

#include "mpi.h"
int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int nprocs;
    int myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();

    return result;
}
#endif