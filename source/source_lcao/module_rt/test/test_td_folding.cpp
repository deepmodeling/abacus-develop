#include "source_base/constants.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_rt/td_folding.h"

#include <gtest/gtest.h>

#include <array>
#include <complex>
#include <map>
#include <vector>

namespace
{
using Complex = std::complex<double>;
using Vector = ModuleBase::Vector3<double>;

// Real and complex HR use distinct instantiations in td_folding.cpp.
template <typename T> T real_space_value(const int row, const int column, const int translation);

template <> double real_space_value<double>(const int row, const int column, const int translation)
{
    return 1.0 + 0.2 * row + 0.3 * column + 0.4 * translation;
}

template <> Complex real_space_value<Complex>(const int row, const int column, const int translation)
{
    const double real = real_space_value<double>(row, column, translation);
    const double imaginary = 0.1 * (row - column) + 0.2 * translation;
    return Complex(real, imaginary);
}

class TdFoldingTest : public testing::Test
{
  protected:
    template <typename T> void check_folding(const Vector &field, const int layout)
    {
        Parallel_Orbitals orbitals;
        orbitals.para_2d().set_serial(4, 4);
        const std::array<int, 2> atom_start = {{0, 2}};
        orbitals.set_atomic_trace(atom_start.data(), 2, 4);
        const std::array<ModuleBase::Vector3<int>, 3> translations = {{{0, 0, 0}, {1, -1, 2}, {-2, 1, -1}}};
        hamilt::HContainer<T> real_space(&orbitals);
        for (int row_atom = 0; row_atom < 2; ++row_atom)
        {
            for (int column_atom = 0; column_atom < 2; ++column_atom)
            {
                for (const auto &translation : translations)
                {
                    hamilt::AtomPair<T> pair(
                        row_atom,
                        column_atom,
                        translation,
                        &orbitals);
                    real_space.insert_pair(pair);
                }
            }
        }
        real_space.allocate(nullptr, true);
        for (int index = 0; index < real_space.size_atom_pairs(); ++index)
        {
            auto &pair = real_space.get_atom_pair(index);
            for (int ir = 0; ir < pair.get_R_size(); ++ir)
            {
                for (int row = 0; row < pair.get_row_size(); ++row)
                {
                    for (int column = 0; column < pair.get_col_size(); ++column)
                    {
                        const int global_row = pair.get_begin_row() + row;
                        const int global_column = pair.get_begin_col() + column;
                        const int element = row * pair.get_col_size() + column;
                        pair.get_pointer(ir)[element] = real_space_value<T>(global_row, global_column, ir);
                    }
                }
            }
        }
        // Explicit Cartesian components provide an oracle for a skew cell, lat0=3.
        std::map<ModuleBase::Vector3<int>, Complex> phases;
        std::array<double, 3> field_angles;
        for (int ir = 0; ir < 3; ++ir)
        {
            const auto &r = translations[ir];
            const double x = 3.0 * r.x;
            const double y = 3.0 * (0.2 * r.x + r.y);
            const double z = 3.0 * (0.1 * r.x + 0.3 * r.y + r.z);
            field_angles[ir] = field.x * x + field.y * y + field.z * z;
            const Complex exponent(0.0, field_angles[ir]);
            phases[r] = std::exp(exponent);
        }
        const std::array<Vector, 2> points = {{{0.0, 0.0, 0.0}, {0.25, -0.125, 0.375}}};
        for (const Vector &kpoint : points)
        {
            // Padding and a nonzero seed also check leading dimension and accumulation.
            const int leading_dimension = 6;
            const Complex seed(0.3, -0.2);
            std::vector<Complex> folded(24, seed);
            module_rt::folding_HR_td(
                real_space,
                folded.data(),
                kpoint,
                field,
                phases,
                leading_dimension,
                layout);
            std::vector<Complex> expected(24, seed);
            for (int row = 0; row < 4; ++row)
            {
                for (int column = 0; column < 4; ++column)
                {
                    const int element =
                        layout == 0 ? row * leading_dimension + column : column * leading_dimension + row;
                    for (int ir = 0; ir < 3; ++ir)
                    {
                        const auto &r = translations[ir];
                        const double bloch_angle =
                            ModuleBase::TWO_PI * (kpoint.x * r.x + kpoint.y * r.y + kpoint.z * r.z);
                        const Complex exponent(0.0, bloch_angle + field_angles[ir]);
                        expected[element] += real_space_value<T>(row, column, ir) * std::exp(exponent);
                    }
                }
            }
            for (int element = 0; element < 24; ++element)
            {
                SCOPED_TRACE(element);
                const double error = std::abs(folded[element] - expected[element]);
                EXPECT_LT(error, 1.0e-12);
            }
        }
    }
};

TEST_F(TdFoldingTest, RealBlocksMatchAnalyticBlochAndHybridPhases)
{
    const Vector no_field;
    const Vector saved_field(0.07, -0.11, 0.05);
    for (const int layout : {0, 1})
    {
        check_folding<double>(no_field, layout);
        check_folding<double>(saved_field, layout);
    }
}

TEST_F(TdFoldingTest, ComplexBlocksMatchAnalyticBlochAndHybridPhases)
{
    const Vector no_field;
    const Vector saved_field(0.07, -0.11, 0.05);
    for (const int layout : {0, 1})
    {
        check_folding<Complex>(no_field, layout);
        check_folding<Complex>(saved_field, layout);
    }
}
} // namespace
