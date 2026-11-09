#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "csr_test_helpers.h"

#include "source_io/module_hs/pos_op_csr.h"

#include <fstream>
#include <string>
#include <vector>

TEST(PosOpCsr, LatRNonemptySkipsEmptyBlocks)
{
    int empty_counts[3] = {0, 0, 0};
    int x_only_counts[3] = {1, 0, 0};
    int z_only_counts[3] = {0, 0, 2};

    EXPECT_FALSE(ModuleIO::detail::lat_r_nonempty(empty_counts));
    EXPECT_TRUE(ModuleIO::detail::lat_r_nonempty(x_only_counts));
    EXPECT_TRUE(ModuleIO::detail::lat_r_nonempty(z_only_counts));
}

TEST(PosOpCsr, CsrAssembleTextAllowsZeroBlocks)
{
    const std::string payload_filename = "rr_empty_payload.tmp";
    const std::string output_filename = "rr_empty.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str());
    payload.close();

    UnitCell ucell;
    init_unitcell(ucell);

    ModuleIO::detail::assemble_csr(output_filename,
                                   payload_filename,
                                   9,
                                   2,
                                   0,
                                   false,
                                   false,
                                   "x",
                                   &ucell,
                                   "PosOpCsr");

    const std::string output = read_file(output_filename);
    EXPECT_THAT(output, testing::HasSubstr("Ionic_Step 9\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print r matrix in real space r(R) x-component\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of spin directions\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # spin index\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 2 # number of localized basis\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 0 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr("user_defined_lattice\n"));
    EXPECT_THAT(output, testing::HasSubstr("Direct\n 0 0.25 0.5\n"));
    EXPECT_THAT(output, testing::HasSubstr("CSR Format"));

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

TEST(PosOpCsr, CsrAssembleTextKeepsSingleDirectionPayload)
{
    const std::string payload_filename = "rr_single_direction_payload.tmp";
    const std::string output_filename = "rr_single_direction.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str());
    payload << "1 0 -1\n";
    payload << "1\n";
    payload << " 4.00000000e+00\n";
    payload << " 0\n";
    payload << "0 1\n";
    payload << "0\n";
    payload << "0\n";
    payload.close();

    UnitCell ucell;
    init_unitcell(ucell);

    ModuleIO::detail::assemble_csr(output_filename,
                                   payload_filename,
                                   10,
                                   2,
                                   1,
                                   false,
                                   false,
                                   "y",
                                   &ucell,
                                   "PosOpCsr");

    const std::string output = read_file(output_filename);
    EXPECT_THAT(output, testing::HasSubstr("Ionic_Step 10\n"));
    EXPECT_THAT(output, testing::HasSubstr(" # print r matrix in real space r(R) y-component\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 1 # number of Bravais lattice vector R\n"));
    EXPECT_THAT(output, testing::HasSubstr("1 0 -1\n"));
    EXPECT_THAT(output, testing::HasSubstr(" 4.00000000e+00\n"));

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

TEST(PosOpCsr, CsrAssembleBinaryKeepsHeaderAndPayloadOrder)
{
    const std::string payload_filename = "rr_binary_payload.tmp";
    const std::string output_filename = "rr_binary.csr";
    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());

    std::ofstream payload(payload_filename.c_str(), std::ios::binary);
    int dRx = 1;
    int dRy = 2;
    int dRz = 3;
    int x_count = 1;
    int y_count = 0;
    int z_count = 0;
    double value = 4.0;
    int column = 1;
    long long ptr0 = 0;
    long long ptr1 = 1;
    payload.write(reinterpret_cast<const char*>(&dRx), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&dRy), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&dRz), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&x_count), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&value), sizeof(double));
    payload.write(reinterpret_cast<const char*>(&column), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&ptr0), sizeof(long long));
    payload.write(reinterpret_cast<const char*>(&ptr1), sizeof(long long));
    payload.write(reinterpret_cast<const char*>(&y_count), sizeof(int));
    payload.write(reinterpret_cast<const char*>(&z_count), sizeof(int));
    payload.close();

    ModuleIO::detail::assemble_csr(output_filename,
                                   payload_filename,
                                   11,
                                   2,
                                   1,
                                   true,
                                   false,
                                   "x",
                                   nullptr,
                                   "PosOpCsr");

    std::ifstream ifs(output_filename, std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    EXPECT_EQ(read_binary_value<int>(ifs), 11);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 2);
    EXPECT_EQ(read_binary_value<int>(ifs), 3);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_DOUBLE_EQ(read_binary_value<double>(ifs), 4.0);
    EXPECT_EQ(read_binary_value<int>(ifs), 1);
    EXPECT_EQ(read_binary_value<long long>(ifs), 0);
    EXPECT_EQ(read_binary_value<long long>(ifs), 1);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    EXPECT_EQ(read_binary_value<int>(ifs), 0);
    ifs.close();

    std::remove(payload_filename.c_str());
    std::remove(output_filename.c_str());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
