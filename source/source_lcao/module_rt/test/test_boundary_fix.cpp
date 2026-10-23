#include "source_lcao/module_rt/boundary_fix.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_cell/mdcell.h"
#include "source_cell/update_cell.h"
#include "source_md/md_base.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

namespace
{
using Complex = std::complex<double>;

// Two artificial atoms, one orbital each, and Gamma plus a non-Gamma point.
// No structure, pseudopotential, orbitals, or SCF reference data are needed.
class BoundaryFixTest : public testing::Test
{
  protected:
    UnitCell cell;
    std::vector<Atom> species;
    K_Vectors kpoints;
    Parallel_Orbitals orbitals;
    const std::array<Complex, 4> initial_h = {{2.0, Complex(0.3, 0.4), Complex(0.3, -0.4), 3.0}};
    const std::array<Complex, 4> initial_s = {{1.0, Complex(0.2, 0.1), Complex(0.2, -0.1), 1.5}};
    const std::array<Complex, 2> initial_c = {{1.0, Complex(0.0, 0.5)}};

    void SetUp() override
    {
        species.resize(1);
        cell.atoms = species.data();
        cell.ntype = 1;
        cell.nat = 2;
        cell.lat0 = 1.0;
        cell.omega = 1.0;
        cell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0);
        cell.GT = cell.latvec;
        Atom& atom = species[0];
        atom.na = 2;
        atom.nw = 1;
        atom.taud.resize(2);
        atom.tau.resize(2);
        atom.dis.resize(2);
        atom.vel.resize(2);
        atom.mbl.resize(2);
        atom.boundary_shift.resize(2);
        atom.taud[0].set(0.0, 0.0, 0.0);
        atom.taud[1].set(0.5, 0.5, 0.5);
        atom.tau = atom.taud;
        cell.set_iat2itia();
        cell.set_iat2iwt(1);
        orbitals.set_serial(2, 2);
        orbitals.ncol_bands = 1;
        const int* atom_offsets = cell.get_iat2iwt();
        orbitals.set_atomic_trace(atom_offsets, 2, 2);
        kpoints.set_nks(2);
        kpoints.kvec_d.resize(2);
        kpoints.kvec_d[0].set(0.0, 0.0, 0.0);
        kpoints.kvec_d[1].set(0.25, 0.125, -0.25);
    }

    Complex quadratic_form(const Complex* matrix, const Complex* coefficients) const
    {
        Complex result = 0.0;
        for (int column = 0; column < 2; ++column)
        {
            for (int row = 0; row < 2; ++row)
            {
                result += std::conj(coefficients[row]) * matrix[row + 2 * column] * coefficients[column];
            }
        }
        return result;
    }

    void expect_phase_correction(const std::array<Complex, 2>& phases)
    {
        const ct::TensorShape matrix_shape({8});
        ct::Tensor h(ct::DataType::DT_COMPLEX_DOUBLE, matrix_shape);
        ct::Tensor s(ct::DataType::DT_COMPLEX_DOUBLE, matrix_shape);
        psi::Psi<Complex> coefficients(2, 1, 2, 2, true);
        for (int ik = 0; ik < 2; ++ik)
        {
            Complex* h_block = h.data<Complex>() + 4 * ik;
            Complex* s_block = s.data<Complex>() + 4 * ik;
            std::copy(initial_h.begin(), initial_h.end(), h_block);
            std::copy(initial_s.begin(), initial_s.end(), s_block);
            coefficients.fix_k(ik);
            std::copy(initial_c.begin(), initial_c.end(), coefficients.get_pointer());
        }
        module_rt::reset_matrix_boundary(cell, kpoints, &orbitals, h, s, &coefficients, 4);
        for (int ik = 0; ik < 2; ++ik)
        {
            const Complex* h_block = h.data<Complex>() + 4 * ik;
            const Complex* s_block = s.data<Complex>() + 4 * ik;
            coefficients.fix_k(ik);
            const Complex* c_block = coefficients.get_pointer();
            for (int column = 0; column < 2; ++column)
            {
                const Complex column_phase = ik == 0 ? Complex(1.0, 0.0) : phases[column];
                const Complex expected_c = std::conj(column_phase) * initial_c[column];
                const double c_error = std::abs(c_block[column] - expected_c);
                EXPECT_NEAR(c_error, 0.0, 1.0e-13);
                for (int row = 0; row < 2; ++row)
                {
                    const Complex row_phase = ik == 0 ? Complex(1.0, 0.0) : phases[row];
                    const int index = row + 2 * column;
                    const Complex factor = std::conj(row_phase) * column_phase;
                    const Complex expected_h = initial_h[index] * factor;
                    const Complex expected_s = initial_s[index] * factor;
                    const double h_error = std::abs(h_block[index] - expected_h);
                    const double s_error = std::abs(s_block[index] - expected_s);
                    EXPECT_NEAR(h_error, 0.0, 1.0e-13);
                    EXPECT_NEAR(s_error, 0.0, 1.0e-13);
                }
            }
            const Complex original_norm = quadratic_form(initial_s.data(), initial_c.data());
            const Complex original_energy = quadratic_form(initial_h.data(), initial_c.data());
            const Complex corrected_norm = quadratic_form(s_block, c_block);
            const Complex corrected_energy = quadratic_form(h_block, c_block);
            const double norm_error = std::abs(corrected_norm - original_norm);
            const double energy_error = std::abs(corrected_energy - original_energy);
            EXPECT_NEAR(norm_error, 0.0, 1.0e-13);
            EXPECT_NEAR(energy_error, 0.0, 1.0e-13);
        }

        // A subsequent adjustment with no crossing must clear the image shift,
        // so the previous step's phase cannot be applied a second time.
        std::array<Complex, 8> corrected_h;
        std::array<Complex, 8> corrected_s;
        std::array<Complex, 4> corrected_c;
        std::copy_n(h.data<Complex>(), 8, corrected_h.begin());
        std::copy_n(s.data<Complex>(), 8, corrected_s.begin());
        for (int ik = 0; ik < 2; ++ik)
        {
            coefficients.fix_k(ik);
            Complex* destination = corrected_c.data() + 2 * ik;
            std::copy_n(coefficients.get_pointer(), 2, destination);
        }
        unitcell::periodic_boundary_adjustment(cell.atoms, cell.latvec, cell.ntype);
        module_rt::reset_matrix_boundary(cell, kpoints, &orbitals, h, s, &coefficients, 4);
        for (int index = 0; index < 8; ++index)
        {
            EXPECT_EQ(h.data<Complex>()[index], corrected_h[index]);
            EXPECT_EQ(s.data<Complex>()[index], corrected_s[index]);
        }
        for (int ik = 0; ik < 2; ++ik)
        {
            coefficients.fix_k(ik);
            for (int basis = 0; basis < 2; ++basis)
            {
                const int index = 2 * ik + basis;
                EXPECT_EQ(coefficients.get_pointer()[basis], corrected_c[index]);
            }
        }
    }
};

TEST_F(BoundaryFixTest, TinyNegativeMdDisplacementHasNoSpuriousImageOrBlochPhase)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    std::vector<LocalAtom> owned;
    if (domain.rank() == 0)
    {
        owned.resize(2);
        for (int index = 0; index < 2; ++index)
        {
            owned[index].type_index = index;
            owned[index].frac = species[0].taud[index];
            owned[index].cart = owned[index].frac;
        }
    }
    const std::vector<std::string> labels(1, "X");
    const std::vector<double> masses(1, 1.0);
    const std::vector<std::int64_t> counts(1, 2);
    MDCell moving_cell;
    moving_cell.initialize_from_owned_atoms(cell.latvec, cell.GT, 1.0, 1.0, 2,
                                           owned, labels, masses, counts, 0.0, domain);
    moving_cell.set_backing_unitcell(cell);
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = 0.0;
    parameters.md_tlast = 0.0;
    MD_base integrator(parameters, false, true, 0, moving_cell);
    if (!moving_cell.owned_atoms().empty())
    {
        moving_cell.owned_atoms()[0].vel.set(-1.0e-21, -2.0e-21, -3.0e-21);
    }
    std::ofstream log;
    integrator.first_half(log);
    moving_cell.sync_backing_unitcell();
    // Check before another boundary adjustment could hide an invalid 1.0.
    for (int axis = 0; axis < 3; ++axis)
    {
        EXPECT_DOUBLE_EQ(species[0].taud[0][axis], 0.0);
        EXPECT_EQ(species[0].boundary_shift[0][axis], 0);
    }
    const std::array<Complex, 2> identity_phases = {{1.0, 1.0}};
    expect_phase_correction(identity_phases);
}

TEST_F(BoundaryFixTest, BothCrossingDirectionsCorrectMatricesAndCoefficientsForEachAxis)
{
    // exp(-2*pi*i*k_axis) at k=(1/4,1/8,-1/4), independent of production sincos.
    const double half_sqrt_two = std::sqrt(0.5);
    const std::array<Complex, 3> positive_shift_phases = {{Complex(0.0, -1.0),
                                                       Complex(half_sqrt_two, -half_sqrt_two),
                                                       Complex(0.0, 1.0)}};
    for (int atom = 0; atom < 2; ++atom)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            for (int shift = -1; shift <= 1; shift += 2)
            {
                SCOPED_TRACE(atom);
                SCOPED_TRACE(axis);
                SCOPED_TRACE(shift);
                species[0].taud[0].set(0.5, 0.5, 0.5);
                species[0].taud[1].set(0.5, 0.5, 0.5);
                species[0].taud[atom][axis] = shift == 1 ? -0.125 : 1.125;
                unitcell::periodic_boundary_adjustment(cell.atoms, cell.latvec, cell.ntype);
                EXPECT_EQ(species[0].boundary_shift[atom][axis], shift);
                const double expected_coordinate = shift == 1 ? 0.875 : 0.125;
                EXPECT_DOUBLE_EQ(species[0].taud[atom][axis], expected_coordinate);
                std::array<Complex, 2> phases = {{1.0, 1.0}};
                phases[atom] = shift == 1 ? positive_shift_phases[axis] : std::conj(positive_shift_phases[axis]);
                expect_phase_correction(phases);
            }
        }
    }
}
} // namespace
