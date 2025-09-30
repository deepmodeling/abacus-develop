#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <memory>
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#ifdef __MPI
#include <mpi.h>
#endif // __MPI

#define private public
#include "source_cell/sep.h"
#undef private

class ReadSepTest : public testing::Test
{
  protected:
    std::string output;
    std::unique_ptr<SepPot> read_sep{new SepPot};

    void SetUp() override
    {
        // Initialization default check
        EXPECT_FALSE(read_sep->is_enable);
        EXPECT_DOUBLE_EQ(read_sep->r_in, 0.0);
        EXPECT_DOUBLE_EQ(read_sep->r_out, 0.0);
        EXPECT_DOUBLE_EQ(read_sep->r_power, 20.0);
        EXPECT_DOUBLE_EQ(read_sep->enhence_a, 1.0);
        EXPECT_EQ(read_sep->mesh, 0);
        EXPECT_EQ(read_sep->strip_elec, 0);
        EXPECT_EQ(read_sep->r, nullptr);
        EXPECT_EQ(read_sep->rv, nullptr);
    }

    void TearDown() override
    {
        // Cleaning is done automatically in the destructor
    }
};

TEST_F(ReadSepTest, ReadSep)
{
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
    {
#endif // !__MPI
        std::ifstream ifs;
        ifs.open("./support/F_pbe_50.sep");
        ASSERT_TRUE(ifs.is_open());
        read_sep->read_sep(ifs);
        ifs.close();
        EXPECT_EQ(read_sep->label, "F");
        EXPECT_EQ(read_sep->mesh, 1038);
        EXPECT_EQ(read_sep->xc_type, "pbe");
        EXPECT_EQ(read_sep->strip_elec, 50);

        EXPECT_EQ(read_sep->r[0], 3.4643182373e-06);
        EXPECT_NE(read_sep->r, nullptr);
        EXPECT_NE(read_sep->rv, nullptr);
#ifdef __MPI
    }
#endif // __MPI
}

TEST_F(ReadSepTest, PrintSep)
{
#ifdef __MPI
    if (GlobalV::MY_RANK == 0)
    {
#endif
        // 设置测试数据
        read_sep->label = "F";
        read_sep->xc_type = "pbe";
        read_sep->orbital = "p";
        read_sep->strip_elec = 50;
        read_sep->mesh = 2;
        read_sep->r = new double[2]{0.1, 0.2};
        read_sep->rv = new double[2]{1.0, 2.0};

        // 测试打印功能
        std::ofstream ofs("test_sep.out");
        read_sep->print_sep_info(ofs);
        read_sep->print_sep_vsep(ofs);
        ofs.close();

        // 验证输出文件
        std::ifstream ifs("test_sep.out");
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(ifs, line))
        {
            lines.push_back(line);
        }
        ifs.close();

        EXPECT_THAT(lines, testing::Contains(" label         F"));
        EXPECT_THAT(lines, testing::Contains(" xc            pbe"));
        EXPECT_THAT(lines, testing::Contains(" orbital       p"));
        EXPECT_THAT(lines, testing::Contains(" strip electron50"));
        EXPECT_THAT(lines, testing::Contains(" mesh  2"));

        std::remove("test_sep.out");
#ifdef __MPI
    }
#endif
}

#ifdef __MPI
TEST_F(ReadSepTest, BcastSep)
{
    if (GlobalV::MY_RANK == 0)
    {
        std::ifstream ifs;
        ifs.open("./support/F_pbe_50.sep");
        ASSERT_TRUE(ifs.is_open());
        read_sep->read_sep(ifs);
        ifs.close();
    }
    read_sep->bcast_sep();
    if (GlobalV::MY_RANK != 0)
    {
        EXPECT_EQ(read_sep->label, "F");
        EXPECT_EQ(read_sep->mesh, 1038);
        EXPECT_EQ(read_sep->xc_type, "pbe");
        EXPECT_EQ(read_sep->strip_elec, 50);
        EXPECT_DOUBLE_EQ(read_sep->r[0], 3.4643182373e-06);
        EXPECT_NE(read_sep->r, nullptr);
        EXPECT_NE(read_sep->rv, nullptr);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);

    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
#endif // __MPI
