#ifndef ABACUS_RT_BOUNDARY_TEST_H
#define ABACUS_RT_BOUNDARY_TEST_H

#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_lcao/module_rt/boundary_fix.h"

#include <array>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

namespace rt_test
{
using Complex = std::complex<double>;
using Vector = ModuleBase::Vector3<double>;

inline Complex matrix_value(int row, int column)
{
    const double real = 1.0 + row + column;
    const double imaginary = 0.1 * (row - column);
    return Complex(real, imaginary);
}

inline Complex boundary_phase(
    const UnitCell &cell,
    const Vector &point,
    const Vector &field,
    int atom)
{
    const auto &shift = cell.atoms[0].boundary_shift[atom];
    const Vector direct(shift.x, shift.y, shift.z);
    const Vector cartesian = direct * cell.latvec * cell.lat0;
    const double angle = -ModuleBase::TWO_PI * (point * direct) - field * cartesian;
    const Complex exponent(0.0, angle);
    return std::exp(exponent);
}

// An analytic oracle for boundary_fix; no MD or real-space folding is involved.
inline std::vector<Complex> analytic_matrix(
    const UnitCell &cell,
    const Parallel_Orbitals &layout,
    const Vector &point,
    const Vector &field,
    bool wrapped)
{
    const int rows = layout.para_2d().nrloc();
    const int columns = layout.para_2d().ncloc();
    std::vector<Complex> matrix(rows * columns);
    for (int column = 0; column < columns; ++column)
    {
        const int j = layout.para_2d().icglb(column);
        for (int row = 0; row < rows; ++row)
        {
            const int i = layout.para_2d().irglb(row);
            Complex phase(1.0, 0.0);
            if (wrapped)
            {
                phase = std::conj(boundary_phase(
                            cell,
                            point,
                            field,
                            i)) *
                        boundary_phase(
                            cell,
                            point,
                            field,
                            j);
            }
            matrix[column * rows + row] = matrix_value(i, j) * phase;
        }
    }
    return matrix;
}

// Shared state/assertions for the direct RT units and the MD -> RT contract.
class BoundaryTest : public testing::TestWithParam<bool>
{
  protected:
    ModuleBase::CommunicationDomain domain;
    std::array<Atom, 1> types;
    UnitCell cell;

    void SetUp() override
    {
        domain = ModuleBase::world_comm_domain();
        types[0].na = 4;
        types[0].nw = 1;
        types[0].boundary_shift.resize(4);
        cell.ntype = 1;
        cell.nat = 4;
        cell.atoms = types.data();
        cell.set_iat2itia();
        cell.set_iat2iwt(1);
        cell.latvec = ModuleBase::Matrix3(
            1.0,
            0.2,
            0.1,
            0.0,
            1.0,
            0.3,
            0.0,
            0.0,
            1.0);
        cell.lat0 = 3.0;
    }

    template <typename MatrixBuilder> void check_correction(const Vector &saved_field, MatrixBuilder build_matrix)
    {
        const bool gathered = GetParam();
        Parallel_Orbitals orbitals;
        const auto communicator = domain.communicator();
        const int init_result = orbitals.para_2d().init(
            4,
            4,
            1,
            communicator,
            false);
        ASSERT_EQ(init_result, 0);
        std::ofstream log;
        const int band_result = orbitals.set_nloc_wfc_Eij(2, log, log);
        ASSERT_EQ(band_result, 0);
        orbitals.set_atomic_trace(cell.get_iat2iwt(), cell.nat, 4);
        Parallel_Orbitals full;
        full.para_2d().set_serial(4, 4);
        full.set_atomic_trace(cell.get_iat2iwt(), cell.nat, 4);
        const Parallel_Orbitals &layout = gathered ? full : orbitals;
        const int matrix_size = layout.para_2d().size();
        const bool gathered_branch = matrix_size != orbitals.para_2d().size();
        const bool expected_gathered_branch = gathered && domain.size() > 1;
        ASSERT_EQ(gathered_branch, expected_gathered_branch);
        const ct::TensorShape shape({4, matrix_size});
        ct::Tensor h(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        ct::Tensor s(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        const std::vector<int> plane_wave_counts;
        const int local_rows = orbitals.para_2d().nrloc();
        psi::Psi<Complex> coefficients(
            4,
            orbitals.ncol_bands,
            local_rows,
            plane_wave_counts,
            true);
        const std::vector<Vector> points = {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {-0.25, 0.0, 0.0}, {0.375, 0.0, 0.0}};
        const Vector next_field(-0.13, 0.17, 0.09);
        std::array<std::vector<Complex>, 4> before;
        std::array<std::vector<Complex>, 4> expected;
        for (int k = 0; k < 4; ++k)
        {
            before[k] = build_matrix(
                cell,
                layout,
                points[k],
                saved_field,
                false);
            expected[k] = build_matrix(
                cell,
                layout,
                points[k],
                saved_field,
                true);
            for (int element = 0; element < matrix_size; ++element)
            {
                h.data<Complex>()[k * matrix_size + element] = before[k][element];
                s.data<Complex>()[k * matrix_size + element] = 0.1 * before[k][element];
            }
            if (saved_field.norm2() > 0.0)
            {
                const auto wrong_time = build_matrix(
                    cell,
                    layout,
                    points[k],
                    next_field,
                    true);
                double difference = 0.0;
                for (int element = 0; element < matrix_size; ++element)
                {
                    difference += std::abs(expected[k][element] - wrong_time[element]);
                }
                domain.sum(&difference, 1);
                EXPECT_GT(difference, 0.01);
            }
            const double real = 1.0 + k;
            const double imaginary = 0.5 - k;
            const Complex initial(real, imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < local_rows; ++row)
                {
                    coefficients(k, column, row) = initial;
                }
            }
        }
        coefficients.fix_k(2);
        module_rt::reset_matrix_boundary(
            cell,
            points,
            4,
            saved_field,
            &orbitals,
            h,
            s,
            &coefficients,
            matrix_size);
        EXPECT_EQ(coefficients.get_current_k(), 2);
        for (int k = 0; k < 4; ++k)
        {
            SCOPED_TRACE(k);
            for (int element = 0; element < matrix_size; ++element)
            {
                const bool corrected = !gathered || domain.rank() == 0;
                const Complex target = corrected ? expected[k][element] : before[k][element];
                const int offset = k * matrix_size + element;
                const double h_error = std::abs(h.data<Complex>()[offset] - target);
                const double s_error = std::abs(s.data<Complex>()[offset] - 0.1 * target);
                EXPECT_LT(h_error, 1.0e-12);
                EXPECT_LT(s_error, 1.0e-12);
            }
            const double real = 1.0 + k;
            const double imaginary = 0.5 - k;
            const Complex initial(real, imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < local_rows; ++row)
                {
                    const int global_row = orbitals.para_2d().irglb(row);
                    const Complex target = initial * std::conj(boundary_phase(
                                                         cell,
                                                         points[k],
                                                         saved_field,
                                                         global_row));
                    const double error = std::abs(coefficients(k, column, row) - target);
                    EXPECT_LT(error, 1.0e-12);
                }
            }
        }
    }
};
} // namespace rt_test
#endif
