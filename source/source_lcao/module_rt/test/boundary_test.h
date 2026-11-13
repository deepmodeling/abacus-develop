#ifndef ABACUS_RT_BOUNDARY_TEST_H
#define ABACUS_RT_BOUNDARY_TEST_H

#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_lcao/module_rt/boundary_fix.h"

#include <algorithm>
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

// The matrix/psi setup and checks were moved from test_md_base.cpp.
class BoundaryTest : public testing::TestWithParam<bool>
{
  protected:
    ModuleBase::CommunicationDomain domain;
    std::array<Atom, 1> types;
    UnitCell cell;
    Parallel_Orbitals orbitals;
    Parallel_Orbitals full;
    ct::Tensor h;
    ct::Tensor s;
    psi::Psi<Complex> coefficients;
    std::vector<Complex> expected;
    const std::vector<Vector> points = {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}, {-0.25, 0.0, 0.0}, {0.375, 0.0, 0.0}};

    const Parallel_Orbitals &matrix_layout() const { return GetParam() ? full : orbitals; }

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
        full.para_2d().set_serial(4, 4);
        full.set_atomic_trace(cell.get_iat2iwt(), cell.nat, 4);
        const int matrix_size = matrix_layout().para_2d().size();
        const bool gathered_branch = matrix_size != orbitals.para_2d().size();
        const bool expected_branch = GetParam() && domain.size() > 1;
        ASSERT_EQ(gathered_branch, expected_branch);
        const ct::TensorShape shape({4, matrix_size});
        h = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        s = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        expected.resize(4 * matrix_size);
        const int local_rows = orbitals.para_2d().nrloc();
        coefficients.resize(4, orbitals.ncol_bands, local_rows);
    }

    Complex phase_for(int atom, int k, const Vector &field) const
    {
        const auto &shift = types[0].boundary_shift[atom];
        const Vector direct(shift.x, shift.y, shift.z);
        const Vector cartesian = direct * cell.latvec * cell.lat0;
        const double angle = -ModuleBase::TWO_PI * (points[k] * direct) - field * cartesian;
        const Complex exponent(0.0, angle);
        return std::exp(exponent);
    }

    void expect_phase_correction(const Vector &saved_A)
    {
        const auto &layout = matrix_layout().para_2d();
        const int rows = layout.nrloc();
        const int matrix_size = layout.size();
        const Vector next_A(-0.13, 0.17, 0.09);
        for (int k = 0; k < 4; ++k)
        {
            double difference = 0.0;
            for (int column = 0; column < layout.ncloc(); ++column)
            {
                const int j = layout.icglb(column);
                for (int row = 0; row < rows; ++row)
                {
                    const int i = layout.irglb(row);
                    const int offset = k * matrix_size + column * rows + row;
                    const auto initial = matrix_value(i, j);
                    h.data<Complex>()[offset] = initial;
                    expected[offset] = initial * std::conj(phase_for(i, k, saved_A)) * phase_for(j, k, saved_A);
                    const auto wrong_time = initial * std::conj(phase_for(i, k, next_A)) * phase_for(j, k, next_A);
                    difference += std::abs(expected[offset] - wrong_time);
                }
            }
            if (saved_A.norm2() > 0.0)
            {
                EXPECT_GT(difference, 0.01);
            }
        }
        expect_correction(saved_A);
    }

    void expect_correction(const Vector &saved_A)
    {
        const int matrix_size = matrix_layout().para_2d().size();
        const int rows = orbitals.para_2d().nrloc();
        const int psi_size = rows * orbitals.ncol_bands;
        const Complex *initial_h = h.data<Complex>();
        const int total_size = 4 * matrix_size;
        const std::vector<Complex> before(initial_h, initial_h + total_size);
        for (int element = 0; element < total_size; ++element)
        {
            s.data<Complex>()[element] = 0.1 * before[element];
        }
        for (int k = 0; k < 4; ++k)
        {
            coefficients.fix_k(k);
            const double real = 1.0 + k;
            const double imaginary = 0.5 - k;
            const Complex initial(real, imaginary);
            std::fill_n(coefficients.get_pointer(), psi_size, initial);
        }
        coefficients.fix_k(0);
        module_rt::reset_matrix_boundary(
            cell,
            points,
            4,
            saved_A,
            &orbitals,
            h,
            s,
            &coefficients,
            matrix_size);
        EXPECT_EQ(coefficients.get_current_k(), 0);
        for (int k = 0; k < 4; ++k)
        {
            for (int element = 0; element < matrix_size; ++element)
            {
                const int offset = k * matrix_size + element;
                const bool corrected = !GetParam() || domain.rank() == 0;
                const Complex target = corrected ? expected[offset] : before[offset];
                const double h_error = std::abs(h.data<Complex>()[offset] - target);
                const double s_error = std::abs(s.data<Complex>()[offset] - 0.1 * target);
                EXPECT_LT(h_error, 1.0e-12);
                EXPECT_LT(s_error, 1.0e-12);
            }
            coefficients.fix_k(k);
            const double real = 1.0 + k;
            const double imaginary = 0.5 - k;
            const Complex initial(real, imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < rows; ++row)
                {
                    const int atom = orbitals.para_2d().irglb(row);
                    const Complex target = initial * std::conj(phase_for(atom, k, saved_A));
                    const int element = column * rows + row;
                    const double error = std::abs(coefficients.get_pointer()[element] - target);
                    EXPECT_LT(error, 1.0e-12);
                }
            }
        }
    }
};
} // namespace rt_test
#endif
