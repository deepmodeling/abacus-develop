#include "../rpa_lri_detail.h"

#include "../rpa_lri_bz_io.h"
#include "../rpa_lri_coeff_v1.h"

#include <RI/global/Tensor.h>
#include <array>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using RpaLriDetail::BzSampleData;
using RpaLriDetail::BzSamplePoint;
using RpaLriDetail::checked_add_u64;
using RpaLriDetail::checked_i32_from_int;
using RpaLriDetail::checked_i32_from_size;
using RpaLriDetail::checked_i64_from_u64;
using RpaLriDetail::checked_mul_u64;
using RpaLriDetail::coulomb_atom_pair_index;
using RpaLriDetail::has_valid_matrix_shape;
using RpaLriDetail::LIBRPA_LRICOEF_V1_MARKER;
using RpaLriDetail::LriCoefV1Block;
using RpaLriDetail::real_as_double;
using RpaLriDetail::sum_int_vector;
using RpaLriDetail::write_bz_sample;
using RpaLriDetail::write_lri_coef_v1;
}

TEST(RpaLriDetail, CoulombAtomPairIndexUsesUpperTriangleOrder)
{
    const std::size_t natoms = 4;
    const std::vector<std::size_t> expected = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::size_t index = 0;
    for (std::size_t I = 0; I < natoms; ++I)
    {
        for (std::size_t J = I; J < natoms; ++J)
        {
            EXPECT_EQ(coulomb_atom_pair_index(I, J, natoms), expected[index]);
            ++index;
        }
    }
    EXPECT_THROW(coulomb_atom_pair_index(2, 1, natoms), std::runtime_error);
}

TEST(RpaLriDetail, CheckedUnsignedArithmeticRejectsOverflow)
{
    const unsigned long long max_value = std::numeric_limits<unsigned long long>::max();
    EXPECT_EQ(checked_mul_u64(6, 7, "product"), 42);
    EXPECT_EQ(checked_add_u64(6, 7, "sum"), 13);
    EXPECT_THROW(checked_mul_u64(max_value, 2, "product"), std::runtime_error);
    EXPECT_THROW(checked_add_u64(max_value, 1, "sum"), std::runtime_error);
}

TEST(RpaLriDetail, CheckedIntegerConversionsRejectOverflowAndNegativeValues)
{
    EXPECT_EQ(checked_i64_from_u64(17, "value"), 17);
    EXPECT_EQ(checked_i32_from_size(17, "value"), 17);
    EXPECT_EQ(checked_i32_from_int(17, "value"), 17);
    EXPECT_THROW(checked_i64_from_u64(std::numeric_limits<unsigned long long>::max(), "value"),
                 std::runtime_error);
    EXPECT_THROW(checked_i32_from_size(std::numeric_limits<std::size_t>::max(), "value"), std::runtime_error);
    EXPECT_THROW(checked_i32_from_int(-1, "value"), std::runtime_error);
}

TEST(RpaLriDetail, SumIntVectorRejectsPositiveOverflow)
{
    EXPECT_EQ(sum_int_vector({1, 2, 3}), 6);
    EXPECT_EQ(sum_int_vector({}), 0);
    EXPECT_THROW(sum_int_vector({std::numeric_limits<int>::max(), 1}), std::runtime_error);
}

TEST(RpaLriDetail, RealPartConversionAndMatrixShapeValidation)
{
    EXPECT_DOUBLE_EQ(real_as_double(2.5), 2.5);
    EXPECT_DOUBLE_EQ(real_as_double(std::complex<double>(2.5, -4.0)), 2.5);

    RI::Tensor<double> tensor;
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
    tensor.shape = {2, 3};
    EXPECT_TRUE(has_valid_matrix_shape(tensor));
    tensor.shape = {2};
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
    tensor.shape = {0, 3};
    EXPECT_FALSE(has_valid_matrix_shape(tensor));
}

TEST(RpaLriDetail, CoefficientV1WritesHeaderOffsetsAndPayload)
{
    const std::string filename = "rpa_lri_coeff_v1_test.dat";
    std::remove(filename.c_str());

    RI::Tensor<double> coeff({2, 1, 2});
    coeff(0, 0, 0) = 1.25;
    coeff(1, 0, 0) = -2.0;
    coeff(0, 0, 1) = 3.5;
    coeff(1, 0, 1) = 4.25;

    LriCoefV1Block block;
    block.ia1 = 2;
    block.ia2 = 1;
    block.r[0] = -1;
    block.r[1] = 2;
    block.r[2] = 0;
    block.max_abs = 4.25;
    block.nw1 = 1;
    block.nw2 = 2;
    block.naux = 2;
    block.tensor = &coeff;
    RI::Tensor<double> coeff2({1, 1, 1});
    coeff2(0, 0, 0) = -7.0;
    LriCoefV1Block block2;
    block2.ia1 = 1;
    block2.ia2 = 1;
    block2.max_abs = 7.0;
    block2.nw1 = 1;
    block2.nw2 = 1;
    block2.naux = 1;
    block2.tensor = &coeff2;
    std::vector<LriCoefV1Block> blocks;
    blocks.push_back(block);
    blocks.push_back(block2);

    write_lri_coef_v1(filename, 2, blocks);

    std::ifstream input(filename.c_str(), std::ios::binary);
    ASSERT_TRUE(input.good());
    std::int32_t marker = 0;
    std::int32_t natom = 0;
    std::int32_t ncell = -1;
    std::int64_t nblocks = 0;
    std::int64_t nrecords = 0;
    input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
    input.read(reinterpret_cast<char*>(&natom), sizeof(natom));
    input.read(reinterpret_cast<char*>(&ncell), sizeof(ncell));
    input.read(reinterpret_cast<char*>(&nblocks), sizeof(nblocks));
    input.read(reinterpret_cast<char*>(&nrecords), sizeof(nrecords));
    ASSERT_TRUE(input.good());
    EXPECT_EQ(marker, LIBRPA_LRICOEF_V1_MARKER);
    EXPECT_EQ(natom, 2);
    EXPECT_EQ(ncell, 0);
    EXPECT_EQ(nblocks, 2);
    EXPECT_EQ(nrecords, 2);

    std::int32_t atom1 = 0;
    std::int32_t atom2 = 0;
    std::int32_t r[3] = {0, 0, 0};
    double max_abs = 0.0;
    std::int64_t offset = 0;
    input.read(reinterpret_cast<char*>(&atom1), sizeof(atom1));
    input.read(reinterpret_cast<char*>(&atom2), sizeof(atom2));
    input.read(reinterpret_cast<char*>(r), sizeof(r));
    input.read(reinterpret_cast<char*>(&max_abs), sizeof(max_abs));
    input.read(reinterpret_cast<char*>(&offset), sizeof(offset));
    ASSERT_TRUE(input.good());
    EXPECT_EQ(atom1, 2);
    EXPECT_EQ(atom2, 1);
    EXPECT_EQ(r[0], -1);
    EXPECT_EQ(r[1], 2);
    EXPECT_EQ(r[2], 0);
    EXPECT_DOUBLE_EQ(max_abs, 4.25);
    EXPECT_EQ(offset, 100);

    std::int32_t atom3 = 0;
    std::int32_t atom4 = 0;
    std::int32_t r2[3] = {0, 0, 0};
    double max_abs2 = 0.0;
    std::int64_t offset2 = 0;
    input.read(reinterpret_cast<char*>(&atom3), sizeof(atom3));
    input.read(reinterpret_cast<char*>(&atom4), sizeof(atom4));
    input.read(reinterpret_cast<char*>(r2), sizeof(r2));
    input.read(reinterpret_cast<char*>(&max_abs2), sizeof(max_abs2));
    input.read(reinterpret_cast<char*>(&offset2), sizeof(offset2));
    ASSERT_TRUE(input.good());
    EXPECT_EQ(atom3, 1);
    EXPECT_EQ(atom4, 1);
    EXPECT_EQ(r2[0], 0);
    EXPECT_EQ(r2[1], 0);
    EXPECT_EQ(r2[2], 0);
    EXPECT_DOUBLE_EQ(max_abs2, 7.0);
    EXPECT_EQ(offset2, 132);

    double values[4] = {};
    input.read(reinterpret_cast<char*>(values), sizeof(values));
    ASSERT_TRUE(input.good());
    EXPECT_DOUBLE_EQ(values[0], 1.25);
    EXPECT_DOUBLE_EQ(values[1], -2.0);
    EXPECT_DOUBLE_EQ(values[2], 3.5);
    EXPECT_DOUBLE_EQ(values[3], 4.25);
    double value2 = 0.0;
    input.read(reinterpret_cast<char*>(&value2), sizeof(value2));
    ASSERT_TRUE(input.good());
    EXPECT_DOUBLE_EQ(value2, -7.0);
    EXPECT_EQ(input.peek(), std::char_traits<char>::eof());
    input.close();
    std::remove(filename.c_str());
}

TEST(RpaLriDetail, CoefficientV1RejectsShapeMismatch)
{
    RI::Tensor<double> coeff({1, 1, 1});
    LriCoefV1Block block;
    block.nw1 = 2;
    block.nw2 = 1;
    block.naux = 1;
    block.tensor = &coeff;
    std::vector<LriCoefV1Block> blocks(1, block);

    EXPECT_THROW(write_lri_coef_v1("unused.dat", 1, blocks), std::runtime_error);
}

namespace
{
BzSamplePoint make_bz_point(const double weight,
                            const std::array<double, 3>& direct,
                            const std::array<double, 3>& cartesian)
{
    BzSamplePoint point;
    point.weight = weight;
    point.direct = direct;
    point.cartesian = cartesian;
    return point;
}

BzSampleData make_reduced_bz_sample()
{
    BzSampleData data;
    data.mesh = {{2, 1, 2}};
    data.reduced.push_back(make_bz_point(2.0, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}));
    data.reduced.push_back(make_bz_point(2.0, {{0.5, 0.0, 0.5}}, {{1.0, 0.0, 1.0}}));
    data.full_cartesian = {{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}, {{1.0, 0.0, 0.0}}, {{1.0, 0.0, 1.0}}};
    data.full_map = {0, 0, 1, 1};
    return data;
}
} // namespace

TEST(RpaLriDetail, BzSampleWritesReducedGridAndFullMap)
{
    const std::string filename = "rpa_lri_bz_sample_test.txt";
    std::remove(filename.c_str());
    write_bz_sample(filename, make_reduced_bz_sample());

    std::ifstream input(filename.c_str());
    ASSERT_TRUE(input.good());
    std::string line;
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    EXPECT_EQ(line, "2     1     2");
    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    EXPECT_EQ(line, "2       2");

    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream first_point(line);
    int index = 0;
    double weight = 0.0;
    double direct_x = 0.0;
    double direct_y = 0.0;
    double direct_z = 0.0;
    double cart_x = 0.0;
    double cart_y = 0.0;
    double cart_z = 0.0;
    int scf_index = 0;
    int coulomb_index = 0;
    first_point >> index >> weight >> direct_x >> direct_y >> direct_z >> cart_x >> cart_y >> cart_z >> scf_index
        >> coulomb_index;
    EXPECT_EQ(index, 1);
    EXPECT_DOUBLE_EQ(weight, 0.5);
    EXPECT_DOUBLE_EQ(direct_x, 0.0);
    EXPECT_DOUBLE_EQ(direct_y, 0.0);
    EXPECT_DOUBLE_EQ(direct_z, 0.0);
    EXPECT_DOUBLE_EQ(cart_x, 0.0);
    EXPECT_DOUBLE_EQ(cart_y, 0.0);
    EXPECT_DOUBLE_EQ(cart_z, 0.0);
    EXPECT_EQ(scf_index, 1);
    EXPECT_EQ(coulomb_index, 1);

    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    std::istringstream second_point(line);
    second_point >> index >> weight >> direct_x >> direct_y >> direct_z >> cart_x >> cart_y >> cart_z >> scf_index
        >> coulomb_index;
    EXPECT_EQ(index, 2);
    EXPECT_DOUBLE_EQ(weight, 0.5);
    EXPECT_DOUBLE_EQ(direct_x, 0.5);
    EXPECT_DOUBLE_EQ(direct_y, 0.0);
    EXPECT_DOUBLE_EQ(direct_z, 0.5);
    EXPECT_DOUBLE_EQ(cart_x, 1.0);
    EXPECT_DOUBLE_EQ(cart_y, 0.0);
    EXPECT_DOUBLE_EQ(cart_z, 1.0);
    EXPECT_EQ(scf_index, 2);
    EXPECT_EQ(coulomb_index, 2);

    ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
    EXPECT_EQ(line, "full_kmap 4");
    for (int full_index = 1; full_index <= 4; ++full_index)
    {
        ASSERT_TRUE(static_cast<bool>(std::getline(input, line)));
        std::istringstream full_point(line);
        int point_index = 0;
        int representative = 0;
        double full_x = 0.0;
        double full_y = 0.0;
        double full_z = 0.0;
        full_point >> point_index;
        full_point >> full_x >> full_y >> full_z;
        full_point >> representative;
        EXPECT_EQ(point_index, full_index);
        EXPECT_DOUBLE_EQ(full_x, full_index == 1 || full_index == 2 ? 0.0 : 1.0);
        EXPECT_DOUBLE_EQ(full_y, 0.0);
        EXPECT_DOUBLE_EQ(full_z, full_index == 1 || full_index == 3 ? 0.0 : 1.0);
        EXPECT_EQ(representative, full_index <= 2 ? 1 : 2);
    }
    EXPECT_FALSE(static_cast<bool>(std::getline(input, line)));
    input.close();
    std::remove(filename.c_str());
}

TEST(RpaLriDetail, BzSampleRejectsInvalidRepresentative)
{
    BzSampleData data = make_reduced_bz_sample();
    data.full_map[3] = 2;
    EXPECT_THROW(write_bz_sample("unused.txt", data), std::runtime_error);
}

TEST(RpaLriDetail, BzSampleRejectsInconsistentFoldedWeight)
{
    BzSampleData data = make_reduced_bz_sample();
    data.reduced[0].weight = 3.0;
    data.reduced[1].weight = 1.0;
    EXPECT_THROW(write_bz_sample("unused.txt", data), std::runtime_error);
}

TEST(RpaLriDetail, BzSampleKeepsFullGridFormatWithoutMapTail)
{
    const std::string filename = "rpa_lri_bz_full_test.txt";
    std::remove(filename.c_str());
    BzSampleData data;
    data.mesh = {{1, 1, 2}};
    data.reduced.push_back(make_bz_point(0.5, {{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}));
    data.reduced.push_back(make_bz_point(0.5, {{0.0, 0.0, 0.5}}, {{0.0, 0.0, 1.0}}));
    data.full_cartesian = {{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}};
    data.full_map = {0, 1};
    write_bz_sample(filename, data);

    std::ifstream input(filename.c_str());
    ASSERT_TRUE(input.good());
    std::string line;
    int line_count = 0;
    while (std::getline(input, line))
    {
        ++line_count;
        EXPECT_NE(line, "full_kmap 2");
    }
    EXPECT_EQ(line_count, 4);
    input.close();
    std::remove(filename.c_str());
}

TEST(RpaLriDetail, BzSampleRejectsGridSizeOutsideIndexRange)
{
    BzSampleData data;
    data.mesh = {{std::numeric_limits<int>::max(), 2, 1}};
    EXPECT_THROW(write_bz_sample("unused.txt", data), std::runtime_error);
}

TEST(RpaLriDetail, BzSampleRejectsOutputOpenFailure)
{
    EXPECT_THROW(write_bz_sample("missing_bz_sample_directory/bz_sample.txt", make_reduced_bz_sample()),
                 std::runtime_error);
}
