#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_lcao/module_rt/boundary_fix.h"

#include <gtest/gtest.h>

#include <array>
#include <complex>
#include <fstream>
#include <vector>

namespace
{
using Complex = std::complex<double>;
using Vector = ModuleBase::Vector3<double>;

class BoundaryFixTest : public testing::TestWithParam<bool>
{
  protected:
    ModuleBase::CommunicationDomain domain;
    std::array<Atom, 1> types;
    UnitCell cell;
    const std::array<double, 4> kx = {{0.0,
        0.25,
        -0.25,
        0.375}};

    void SetUp() override
    {
        domain = ModuleBase::world_comm_domain();
        types[0].na = 4;
        types[0].nw = 1;
        types[0].boundary_shift.resize(4);
        types[0].boundary_shift[0].set(-1, 1, -1);
        types[0].boundary_shift[1].set(1, 0, 1);
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

    Complex phase_for(const int atom, const int kpoint, const Vector &field) const
    {
        const auto &shift = types[0].boundary_shift[atom];
        const Vector direct(shift.x, shift.y, shift.z);
        const Vector cartesian = direct * cell.latvec * cell.lat0;
        const double angle = -ModuleBase::TWO_PI * kx[kpoint] * shift.x - field * cartesian;
        const Complex exponent(0.0, angle);
        return std::exp(exponent);
    }

    Complex matrix_value(const int row, const int column) const
    {
        const double real = 1.0 + row + column;
        const double imaginary = 0.1 * (row - column);
        return Complex(real, imaginary);
    }

    void check_correction(const Vector &saved_field)
    {
        const bool gathered = GetParam();
        const int basis_size = 4;
        const int band_count = 2;
        Parallel_Orbitals orbitals;
        const auto communicator = domain.communicator();
        const int init_result = orbitals.para_2d().init(
            basis_size,
            basis_size,
            1,
            communicator,
            false);
        ASSERT_EQ(init_result, 0);
        std::ofstream log;
        const int band_result = orbitals.set_nloc_wfc_Eij(band_count, log, log);
        ASSERT_EQ(band_result, 0);
        orbitals.set_atomic_trace(cell.get_iat2iwt(), cell.nat, basis_size);
        const int rows = gathered ? basis_size : orbitals.para_2d().nrloc();
        const int columns = gathered ? basis_size : orbitals.para_2d().ncloc();
        const int matrix_size = rows * columns;
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
        std::vector<Vector> points(4);
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            points[kpoint].set(kx[kpoint], 0.0, 0.0);
            for (int column = 0; column < columns; ++column)
            {
                const int global_column = gathered ? column : orbitals.para_2d().icglb(column);
                for (int row = 0; row < rows; ++row)
                {
                    const int global_row = gathered ? row : orbitals.para_2d().irglb(row);
                    const int element = kpoint * matrix_size + column * rows + row;
                    h.data<Complex>()[element] = matrix_value(global_row, global_column);
                    s.data<Complex>()[element] = 0.1 * h.data<Complex>()[element];
                }
            }
            const double real = 1.0 + kpoint;
            const double imaginary = 0.5 - kpoint;
            const Complex initial(real, imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < local_rows; ++row)
                {
                    coefficients(kpoint, column, row) = initial;
                }
            }
        }
        // Use a nonzero cursor to detect an accidental fix_k() in the OpenMP loop.
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
        const Vector next_field(-0.13, 0.17, 0.09);
        double field_time_difference = 0.0;
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            SCOPED_TRACE(kpoint);
            for (int column = 0; column < columns; ++column)
            {
                const int global_column = gathered ? column : orbitals.para_2d().icglb(column);
                for (int row = 0; row < rows; ++row)
                {
                    const int global_row = gathered ? row : orbitals.para_2d().irglb(row);
                    const int element = kpoint * matrix_size + column * rows + row;
                    const Complex initial = matrix_value(global_row, global_column);
                    Complex phase(1.0, 0.0);
                    if (!gathered || domain.rank() == 0)
                    {
                        phase = std::conj(phase_for(global_row, kpoint, saved_field)) *
                                phase_for(global_column, kpoint, saved_field);
                        const Complex wrong_phase = std::conj(phase_for(global_row, kpoint, next_field)) *
                                                    phase_for(global_column, kpoint, next_field);
                        field_time_difference += std::abs(initial * (phase - wrong_phase));
                    }
                    const Complex expected_h = initial * phase;
                    const Complex expected_s = 0.1 * expected_h;
                    const double h_error = std::abs(h.data<Complex>()[element] - expected_h);
                    const double s_error = std::abs(s.data<Complex>()[element] - expected_s);
                    EXPECT_LT(h_error, 1.0e-12);
                    EXPECT_LT(s_error, 1.0e-12);
                }
            }
            const double real = 1.0 + kpoint;
            const double imaginary = 0.5 - kpoint;
            const Complex initial(real, imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < local_rows; ++row)
                {
                    const int global_row = orbitals.para_2d().irglb(row);
                    const Complex expected = initial * std::conj(phase_for(global_row, kpoint, saved_field));
                    const double error = std::abs(coefficients(kpoint, column, row) - expected);
                    EXPECT_LT(error, 1.0e-12);
                }
            }
        }
        if (saved_field.norm2() > 0.0)
        {
            const double global_difference = domain.sum(field_time_difference);
            EXPECT_GT(global_difference, 0.01);
        }
    }
};

TEST_P(BoundaryFixTest, BlochPhaseCovariance)
{
    const Vector no_field;
    check_correction(no_field);
}

TEST_P(BoundaryFixTest, HybridUsesSavedFieldInNonorthogonalCell)
{
    const Vector saved_field(0.07, -0.11, 0.05);
    check_correction(saved_field);
}

TEST_P(BoundaryFixTest, NoBoundaryShiftPreservesSavedState)
{
    for (auto &shift : types[0].boundary_shift)
    {
        shift.set(0, 0, 0);
    }
    const Vector no_field;
    check_correction(no_field);
}

INSTANTIATE_TEST_SUITE_P(MatrixLayouts, BoundaryFixTest, testing::Bool());
} // namespace
