#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_io/module_hs/pos_op_csr.h"

#include <fstream>
#include <string>
#include <vector>

namespace
{
std::vector<std::string> read_lines(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line))
    {
        lines.push_back(line);
    }
    return lines;
}

template <typename T>
T read_binary_value(std::ifstream& ifs)
{
    T value{};
    ifs.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}
} // namespace

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

    ModuleIO::detail::assemble_csr(output_filename,
                                   payload_filename,
                                   9,
                                   2,
                                   0,
                                   false,
                                   false,
                                   "PosOpCsr");

    const std::vector<std::string> lines = read_lines(output_filename);
    ASSERT_EQ(lines.size(), 3);
    EXPECT_EQ(lines[0], "STEP: 9");
    EXPECT_EQ(lines[1], "Matrix Dimension of r(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of r(R): 0");

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

    ModuleIO::detail::assemble_csr(output_filename,
                                   payload_filename,
                                   10,
                                   2,
                                   1,
                                   false,
                                   false,
                                   "PosOpCsr");

    const std::vector<std::string> lines = read_lines(output_filename);
    ASSERT_EQ(lines.size(), 10);
    EXPECT_EQ(lines[0], "STEP: 10");
    EXPECT_EQ(lines[1], "Matrix Dimension of r(R): 2");
    EXPECT_EQ(lines[2], "Matrix number of r(R): 1");
    EXPECT_EQ(lines[3], "1 0 -1");
    EXPECT_EQ(lines[4], "1");
    EXPECT_EQ(lines[5], " 4.00000000e+00");
    EXPECT_EQ(lines[8], "0");
    EXPECT_EQ(lines[9], "0");

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
