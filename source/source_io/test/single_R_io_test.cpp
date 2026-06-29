#include "gtest/gtest.h"
#include "gmock/gmock.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "source_io/module_hs/single_R_io.h"
#include "source_base/global_variable.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
/************************************************
 *  unit test of output_single_R
 ***********************************************/
/**
 * - Tested Functions:
 *   - ModuleIO::output_single_R
 *     - output single R data
 */
Parallel_Orbitals::Parallel_Orbitals()
{
}

Parallel_Orbitals::~Parallel_Orbitals()
{
}

void Parallel_2D::set_serial(const int M_A, const int N_A)
{
    this->nrow = M_A;
    this->ncol = N_A;
    this->is_serial = true;
    this->global2local_row_.resize(M_A);
    this->global2local_row_[0] = 0;
    this->global2local_row_[1] = 1;
    this->global2local_row_[2] = -1;
    this->global2local_row_[3] = 2;
    this->global2local_row_[4] = -1; //Some rows have global2local_row_ < 0
}

int Parallel_2D::get_global_row_size() const
{
    return this->nrow;
}

TEST(ModuleIOTest, OutputSingleR)
{
    // Create temporary output file
    std::stringstream ofs_filename;
    GlobalV::DRANK=0;
    ofs_filename << "test_output_single_R_" << GlobalV::DRANK << ".dat";
    std::ofstream ofs(ofs_filename.str());

    // Define input parameters
    const double sparse_threshold = 1e-8;
    const bool binary = false;
    Parallel_Orbitals pv;
    PARAM.sys.nlocal = 99;
    pv.set_serial(5, 5);
    std::map<size_t, std::map<size_t, double>> XR = {
        {0, {{1, 0.5}, {3, 0.3}}},
        {1, {{0, 0.2}, {2, 0.4}}},
        {3, {{1, 0.1}, {4, 0.7}}}
    };
    ModuleIO::SparseWriteOptions options;
    options.threshold = sparse_threshold;
    options.binary = binary;
    options.reduce = true;
    options.temp_dir = "./";

    // Call function under test
    ModuleIO::output_single_R(ofs, XR, pv, options);

    // Close output file and open it for reading
    ofs.close();
    std::ifstream ifs;
    ifs.open("test_output_single_R_0.dat");
    std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
    std::istringstream content(str);
    std::string value_line;
    std::string column_line;
    std::string indptr_line;
    std::getline(content, value_line);
    std::getline(content, column_line);
    std::getline(content, indptr_line);

    std::istringstream value_stream(value_line);
    std::vector<double> values;
    double value = 0.0;
    while (value_stream >> value)
    {
        values.push_back(value);
    }
    ASSERT_EQ(values.size(), 6);
    EXPECT_DOUBLE_EQ(values[0], 0.5);
    EXPECT_DOUBLE_EQ(values[1], 0.3);
    EXPECT_DOUBLE_EQ(values[2], 0.2);
    EXPECT_DOUBLE_EQ(values[3], 0.4);
    EXPECT_DOUBLE_EQ(values[4], 0.1);
    EXPECT_DOUBLE_EQ(values[5], 0.7);

    std::istringstream column_stream(column_line);
    std::vector<int> columns;
    int column = 0;
    while (column_stream >> column)
    {
        columns.push_back(column);
    }
    EXPECT_THAT(columns, testing::ElementsAre(1, 3, 0, 2, 1, 4));

    std::istringstream indptr_stream(indptr_line);
    std::vector<long long> indptr;
    long long index = 0;
    while (indptr_stream >> index)
    {
        indptr.push_back(index);
    }
    EXPECT_THAT(indptr, testing::ElementsAre(0, 2, 4, 4, 6, 6));
    std::remove("test_output_single_R_0.dat");
}

TEST(ModuleIOTest, OutputSingleRComplexKeepsHighPrecision)
{
    const std::string filename = "test_output_single_R_complex.dat";
    std::remove(filename.c_str());
    GlobalV::DRANK = 0;
    std::ofstream ofs(filename);

    Parallel_Orbitals pv;
    pv.set_serial(5, 5);
    ModuleIO::SparseRBlock<std::complex<double>> XR = {
        {0, {{1, std::complex<double>(1.234567890123456, -2.345678901234567)}}}
    };
    ModuleIO::SparseWriteOptions options;
    options.threshold = 1e-12;
    options.binary = false;
    options.reduce = false;
    options.temp_dir = "./";

    ModuleIO::output_single_R(ofs, XR, pv, options);
    ofs.close();

    std::ifstream ifs(filename);
    const std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("(1.2345678901234560e+00,-2.3456789012345669e+00)"));
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("(1.23457,-2.34568)")));

    std::remove(filename.c_str());
}

TEST(ModuleIOTest, OutputSingleRUsesConfiguredPrecision)
{
    const std::string filename = "test_output_single_R_precision.dat";
    std::remove(filename.c_str());
    GlobalV::DRANK = 0;
    std::ofstream ofs(filename);

    Parallel_Orbitals pv;
    pv.set_serial(5, 5);
    ModuleIO::SparseRBlock<std::complex<double>> XR = {
        {0, {{1, std::complex<double>(1.234567890123456, -2.5)}}}
    };
    ModuleIO::SparseWriteOptions options;
    options.threshold = 1e-12;
    options.binary = false;
    options.precision = 8;
    options.reduce = false;
    options.temp_dir = "./";

    ModuleIO::output_single_R(ofs, XR, pv, options);
    ofs.close();

    std::ifstream ifs(filename);
    const std::string output((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(output, testing::HasSubstr("(1.23456789e+00,-2.50000000e+00)"));
    EXPECT_THAT(output, testing::Not(testing::HasSubstr("1.2345678901234560e+00")));

    std::remove(filename.c_str());
}

void write_out_of_range_sparse_column(const char* filename)
{
    GlobalV::DRANK = 0;
    std::ofstream ofs(filename);
    Parallel_Orbitals pv;
    pv.set_serial(5, 5);
    ModuleIO::SparseRBlock<double> XR;
    XR[0][5] = 1.0;
    ModuleIO::SparseWriteOptions options;
    options.threshold = 1e-12;
    options.binary = false;
    options.reduce = false;
    options.temp_dir = "/tmp/";
    ModuleIO::output_single_R(ofs, XR, pv, options);
}

TEST(ModuleIOTest, OutputSingleRRejectsOutOfRangeColumn)
{
    const char* filename = "/tmp/test_output_single_R_invalid.dat";
    std::remove(filename);
    EXPECT_EXIT(
        write_out_of_range_sparse_column(filename),
        ::testing::ExitedWithCode(1),
        "Sparse column index out of range");
    std::remove(filename);
}

int main(int argc, char **argv)
{

#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD,&GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
