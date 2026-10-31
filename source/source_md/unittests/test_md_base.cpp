#include "source_md/md_base.h"
#include "source_md/md_func.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_global.h"
#include "source_cell/domain_decomposition.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"
#include "source_lcao/module_rt/boundary_fix.h"
#include "source_lcao/module_rt/td_folding.h"

#include <gtest/gtest.h>

#include <array>
#include <complex>
#include <fstream>
#include <map>
#include <vector>

namespace
{
class PositionStateSolver : public ModuleESolver::ESolver
{
public:
    std::vector<bool> observed_position_updates;

    void before_all_runners(BaseCell&, const Input_para&) override {}
    void after_all_runners(BaseCell&) override {}
    void runner(BaseCell& base, const int) override
    {
        const UnitCell& cell = static_cast<const UnitCell&>(base);
        observed_position_updates.push_back(cell.ionic_position_updated);
    }
    double cal_energy() override { return 0.0; }
    void cal_force(BaseCell&, ModuleBase::matrix& force) override { force.zero_out(); }
    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override { stress.zero_out(); }
};

class MdBoundaryTest : public testing::Test
{
protected:
    ModuleBase::CommunicationDomain domain;
    ModuleBase::Matrix3 lattice;
    std::array<Atom, 1> backing_types;
    std::array<Atom, 1> reference_types;
    UnitCell backing;
    MDCell cell;
    DomainDecomposition decomposition;
    MD_para parameters;

    void SetUp() override
    {
        domain = ModuleBase::world_comm_domain();
        lattice.e11 = 1.0;
        lattice.e22 = 1.0;
        lattice.e33 = 1.0;
        for (Atom* type : {backing_types.data(), reference_types.data()})
        {
            type->na = 4;
            type->nw = 1;
            type->tau.resize(4);
            type->taud.resize(4);
            type->dis.resize(4);
            type->vel.resize(4);
            type->mbl.resize(4, ModuleBase::Vector3<int>(1, 1, 1));
            type->boundary_shift.resize(4);
            type->taud[0].set(0.99, 0.25, 0.75);
            type->taud[1].set(0.01, 0.75, 0.25);
            type->taud[2].set(0.0, 0.0, 0.0);
            type->taud[3].set(0.5, 0.5, 0.5);
            type->mbl[3].set(0, 0, 0);
            for (int index = 0; index < 4; ++index)
            {
                type->tau[index] = type->taud[index] * lattice;
            }
        }
        backing.ntype = 1;
        backing.nat = 4;
        backing.atoms = backing_types.data();
        backing.set_iat2itia();
        backing.set_iat2iwt(1);
        std::vector<LocalAtom> owned;
        const int source_rank = domain.size() - 1;
        if (domain.rank() == source_rank)
        {
            for (int index = 0; index < 4; ++index)
            {
                LocalAtom atom;
                atom.frac = backing_types[0].taud[index];
                atom.cart = backing_types[0].tau[index];
                atom.mbl = backing_types[0].mbl[index];
                atom.type_index = index;
                atom.owner_rank = source_rank;
                owned.push_back(atom);
            }
        }
        const ModuleBase::Matrix3 inverse = lattice.Inverse();
        const std::vector<std::string> labels(1, "X");
        const std::vector<double> masses(1, 1.0);
        const std::vector<std::int64_t> counts(1, 4);
        cell.initialize_from_owned_atoms(lattice, inverse, 1.0, 1.0, 4, owned,
                                         labels, masses, counts, 0.0, domain);
        cell.set_backing_unitcell(backing);
        decomposition.init(domain, lattice, 1.0, 0.0, 0.0);
        cell.set_neighbor_cutoff(0.05);
        parameters.md_dt = ModuleBase::AU_to_FS;
        parameters.md_tfirst = -1.0;
        parameters.md_tlast = -1.0;
    }

    ModuleBase::Vector3<double> velocity_for(int index, int step) const
    {
        if (step != 0)
        {
            return ModuleBase::Vector3<double>(0.0, 0.01, 0.0);
        }
        if (index == 0)
        {
            return ModuleBase::Vector3<double>(0.02, 0.6, 0.0);
        }
        if (index == 1)
        {
            return ModuleBase::Vector3<double>(-0.02, -0.6, 0.0);
        }
        return ModuleBase::Vector3<double>(-1.0e-21, -2.0e-21, -3.0e-21);
    }

    void expect_local_state(const LocalAtom& atom)
    {
        const int index = atom.type_index;
        const Atom& reference = reference_types[0];
        SCOPED_TRACE(index);
        for (int axis = 0; axis < 3; ++axis)
        {
            EXPECT_DOUBLE_EQ(atom.frac[axis], reference.taud[index][axis]);
            EXPECT_DOUBLE_EQ(atom.displacement[axis], reference.dis[index][axis]);
            EXPECT_EQ(atom.boundary_shift[axis], reference.boundary_shift[index][axis]);
        }
    }

    void expect_backing_state(const bool positions_updated)
    {
        const Atom& actual = backing_types[0];
        const Atom& reference = reference_types[0];
        EXPECT_EQ(backing.ionic_position_updated, positions_updated);
        for (int index = 0; index < 4; ++index)
        {
            SCOPED_TRACE(index);
            for (int axis = 0; axis < 3; ++axis)
            {
                EXPECT_DOUBLE_EQ(actual.taud[index][axis], reference.taud[index][axis]);
                EXPECT_DOUBLE_EQ(actual.tau[index][axis], reference.tau[index][axis]);
                EXPECT_DOUBLE_EQ(actual.dis[index][axis], reference.dis[index][axis]);
                EXPECT_EQ(actual.boundary_shift[index][axis], reference.boundary_shift[index][axis]);
            }
        }
    }
    using Complex = std::complex<double>;

    Complex phase_for(int index, int kpoint) const
    {
        const std::array<double, 4> kx_values = {{0.0, 0.25, -0.25, 0.375}};
        const double kx = kx_values[kpoint];
        const double argument = -ModuleBase::TWO_PI * kx * reference_types[0].boundary_shift[index].x;
        const Complex exponent(0.0, argument);
        return std::exp(exponent);
    }

    Complex matrix_value(int row, int column) const
    {
        const double real = 1.0 + row + column;
        const double imaginary = 0.1 * (row - column);
        return Complex(real, imaginary);
    }

    std::vector<Complex> fold_hybrid(const Parallel_Orbitals& orbitals,
                                     const ModuleBase::Vector3<double>& kpoint,
                                     const ModuleBase::Vector3<double>& vector_potential,
                                     const bool wrapped) const
    {
        hamilt::HContainer<Complex> real_space(&orbitals);
        for (int row = 0; row < 4; ++row)
        {
            if (orbitals.get_nrow_atom(row) == 0)
            {
                continue;
            }
            for (int column = 0; column < 4; ++column)
            {
                if (orbitals.get_ncol_atom(column) == 0)
                {
                    continue;
                }
                const int translation_x = column - row;
                const int translation_y = row - column;
                ModuleBase::Vector3<int> translation(translation_x, translation_y, 0);
                if (wrapped)
                {
                    translation += backing_types[0].boundary_shift[row] - backing_types[0].boundary_shift[column];
                }
                hamilt::AtomPair<Complex> pair(row, column, translation, &orbitals);
                real_space.insert_pair(pair);
            }
        }
        real_space.allocate(nullptr, true);
        std::map<ModuleBase::Vector3<int>, Complex> phases;
        for (int index = 0; index < real_space.size_atom_pairs(); ++index)
        {
            hamilt::AtomPair<Complex>& pair = real_space.get_atom_pair(index);
            const ModuleBase::Vector3<int> translation = pair.get_R_index(0);
            const int row = pair.get_atom_i();
            const int column = pair.get_atom_j();
            pair.get_pointer(0)[0] = matrix_value(row, column);
            const ModuleBase::Vector3<double> direct(translation.x, translation.y, translation.z);
            const ModuleBase::Vector3<double> cartesian = direct * backing.latvec * backing.lat0;
            const double angle = vector_potential * cartesian;
            const Complex exponent(0.0, angle);
            phases[translation] = std::exp(exponent);
        }
        std::vector<Complex> folded(orbitals.nloc, Complex(0.0, 0.0));
        const int rows = orbitals.get_row_size();
        module_rt::folding_HR_td(real_space, folded.data(), kpoint, vector_potential, phases, rows, 1);
        return folded;
    }

    void expect_hybrid_fold_correction(const bool gathered)
    {
        const int basis_size = 4;
        const int band_count = 2;
        Parallel_Orbitals orbitals;
        const int init_result = orbitals.init(basis_size, basis_size, 1, cell.communicator());
        ASSERT_EQ(init_result, 0);
        std::ofstream log;
        const int band_result = orbitals.set_nloc_wfc_Eij(band_count, log, log);
        ASSERT_EQ(band_result, 0);
        orbitals.set_atomic_trace(backing.get_iat2iwt(), backing.nat, basis_size);
        Parallel_Orbitals full;
        full.set_serial(basis_size, basis_size);
        full.set_atomic_trace(backing.get_iat2iwt(), backing.nat, basis_size);
        const Parallel_Orbitals& matrix_layout = gathered ? full : orbitals;
        const int matrix_size = matrix_layout.nloc;
        const ct::TensorShape shape({4, matrix_size});
        ct::Tensor h(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        ct::Tensor s(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        const std::vector<int> plane_wave_counts;
        psi::Psi<Complex> coefficients(4, orbitals.ncol_bands, orbitals.nrow, plane_wave_counts, true);
        const std::array<double, 4> kx_values = {{0.0, 0.25, -0.25, 0.375}};
        std::vector<ModuleBase::Vector3<double>> points(4);
        // A nonorthogonal cell and lat0 != 1 catch missing Cartesian/unit conversions.
        backing.latvec = ModuleBase::Matrix3(1.0, 0.2, 0.1, 0.0, 1.0, 0.3, 0.0, 0.0, 1.0);
        backing.lat0 = 3.0;
        const ModuleBase::Vector3<double> saved_A(0.07, -0.11, 0.05);
        const ModuleBase::Vector3<double> next_A(-0.13, 0.17, 0.09);
        std::vector<std::vector<Complex>> before(4);
        std::vector<std::vector<Complex>> expected(4);
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            points[kpoint].set(kx_values[kpoint], 0.0, 0.0);
            before[kpoint] = fold_hybrid(matrix_layout, points[kpoint], saved_A, false);
            expected[kpoint] = fold_hybrid(matrix_layout, points[kpoint], saved_A, true);
            const std::vector<Complex> wrong_time = fold_hybrid(matrix_layout, points[kpoint], next_A, true);
            double field_time_difference = 0.0;
            for (int element = 0; element < matrix_size; ++element)
            {
                const int offset = kpoint * matrix_size + element;
                h.data<Complex>()[offset] = before[kpoint][element];
                s.data<Complex>()[offset] = 0.1 * before[kpoint][element];
                field_time_difference += std::abs(expected[kpoint][element] - wrong_time[element]);
            }
            EXPECT_GT(field_time_difference, 0.01);
            coefficients.fix_k(kpoint);
            const double coefficient_real = 1.0 + kpoint;
            const double coefficient_imaginary = 0.5 - kpoint;
            const Complex initial(coefficient_real, coefficient_imaginary);
            for (int element = 0; element < orbitals.nrow * orbitals.ncol_bands; ++element)
            {
                coefficients.get_pointer()[element] = initial;
            }
        }
        coefficients.fix_k(0);
        module_rt::reset_matrix_boundary(backing, points, 4, saved_A,
                                         &orbitals, h, s, &coefficients, matrix_size);
        EXPECT_EQ(coefficients.get_current_k(), 0);
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            for (int element = 0; element < matrix_size; ++element)
            {
                const int offset = kpoint * matrix_size + element;
                const bool corrected = !gathered || domain.rank() == 0;
                const Complex target = corrected ? expected[kpoint][element] : before[kpoint][element];
                const double h_error = std::abs(h.data<Complex>()[offset] - target);
                const double s_error = std::abs(s.data<Complex>()[offset] - 0.1 * target);
                EXPECT_LT(h_error, 1.0e-12);
                EXPECT_LT(s_error, 1.0e-12);
            }
            coefficients.fix_k(kpoint);
            const double coefficient_real = 1.0 + kpoint;
            const double coefficient_imaginary = 0.5 - kpoint;
            const Complex initial(coefficient_real, coefficient_imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < orbitals.nrow; ++row)
                {
                    const int global_row = orbitals.local2global_row(row);
                    const ModuleBase::Vector3<int>& shift = backing_types[0].boundary_shift[global_row];
                    const ModuleBase::Vector3<double> direct(shift.x, shift.y, shift.z);
                    const ModuleBase::Vector3<double> cartesian = direct * backing.latvec * backing.lat0;
                    const double angle = ModuleBase::TWO_PI * (points[kpoint] * direct) + saved_A * cartesian;
                    const Complex exponent(0.0, angle);
                    const Complex target = initial * std::exp(exponent);
                    const int element = column * orbitals.nrow + row;
                    const double error = std::abs(coefficients.get_pointer()[element] - target);
                    EXPECT_LT(error, 1.0e-12);
                }
            }
        }
    }

    void expect_phase_correction(bool gathered)
    {
        const int basis_size = 4;
        const int band_count = 2;
        Parallel_Orbitals orbitals;
        const int block_size = 1;
        const int init_result = orbitals.init(basis_size, basis_size, block_size, cell.communicator());
        ASSERT_EQ(init_result, 0);
        std::ofstream log;
        const int band_result = orbitals.set_nloc_wfc_Eij(band_count, log, log);
        ASSERT_EQ(band_result, 0);
        orbitals.set_atomic_trace(backing.get_iat2iwt(), backing.nat, basis_size);
        std::vector<ModuleBase::Vector3<double>> points(4);
        points[1].set(0.25, 0.0, 0.0);
        points[2].set(-0.25, 0.0, 0.0);
        points[3].set(0.375, 0.0, 0.0);
        const int rows = gathered ? basis_size : orbitals.get_row_size();
        const int columns = gathered ? basis_size : orbitals.get_col_size();
        const int matrix_size = rows * columns;
        const ct::TensorShape shape({4, matrix_size});
        ct::Tensor h(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        ct::Tensor s(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        const std::vector<int> plane_wave_counts;
        psi::Psi<Complex> coefficients(4, orbitals.ncol_bands, orbitals.nrow, plane_wave_counts, true);
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            for (int column = 0; column < columns; ++column)
            {
                const int global_column = gathered ? column : orbitals.local2global_col(column);
                for (int row = 0; row < rows; ++row)
                {
                    const int global_row = gathered ? row : orbitals.local2global_row(row);
                    const int element = kpoint * matrix_size + column * rows + row;
                    h.data<Complex>()[element] = matrix_value(global_row, global_column);
                    s.data<Complex>()[element] = 0.1 * h.data<Complex>()[element];
                }
            }
            coefficients.fix_k(kpoint);
            const double coefficient_real = 1.0 + kpoint;
            const double coefficient_imaginary = 0.5 - kpoint;
            const Complex initial(coefficient_real, coefficient_imaginary);
            for (int element = 0; element < orbitals.nrow * orbitals.ncol_bands; ++element)
            {
                coefficients.get_pointer()[element] = initial;
            }
        }
        coefficients.fix_k(0);
        const ModuleBase::Vector3<double> no_vector_potential;
        module_rt::reset_matrix_boundary(backing, points, 4, no_vector_potential,
                                         &orbitals, h, s, &coefficients, matrix_size);
        // Per-k coefficient access must not mutate a shared Psi cursor.
        EXPECT_EQ(coefficients.get_current_k(), 0);
        for (int kpoint = 0; kpoint < 4; ++kpoint)
        {
            for (int column = 0; column < columns; ++column)
            {
                const int global_column = gathered ? column : orbitals.local2global_col(column);
                for (int row = 0; row < rows; ++row)
                {
                    const int global_row = gathered ? row : orbitals.local2global_row(row);
                    const int element = kpoint * matrix_size + column * rows + row;
                    Complex phase(1.0, 0.0);
                    if (!gathered || domain.rank() == 0)
                    {
                        phase = std::conj(phase_for(global_row, kpoint)) * phase_for(global_column, kpoint);
                    }
                    const Complex expected_h = matrix_value(global_row, global_column) * phase;
                    const Complex expected_s = 0.1 * expected_h;
                    const double h_error = std::abs(h.data<Complex>()[element] - expected_h);
                    const double s_error = std::abs(s.data<Complex>()[element] - expected_s);
                    EXPECT_LT(h_error, 1.0e-12);
                    EXPECT_LT(s_error, 1.0e-12);
                }
            }
            coefficients.fix_k(kpoint);
            const double coefficient_real = 1.0 + kpoint;
            const double coefficient_imaginary = 0.5 - kpoint;
            const Complex initial(coefficient_real, coefficient_imaginary);
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < orbitals.nrow; ++row)
                {
                    const int global_row = orbitals.local2global_row(row);
                    const Complex expected = initial * std::conj(phase_for(global_row, kpoint));
                    const int element = column * orbitals.nrow + row;
                    const double error = std::abs(coefficients.get_pointer()[element] - expected);
                    EXPECT_LT(error, 1.0e-12);
                }
            }
        }
    }
};

TEST_F(MdBoundaryTest, MatchesLegacyWrappingThroughMigrationAndRepeatedSync)
{
    const int rank = domain.rank();
    MD_base integrator(parameters, false, true, rank, cell);
    std::ofstream log;
    for (int step = 0; step < 2; ++step)
    {
        SCOPED_TRACE(step);
        for (int index = 0; index < 4; ++index)
        {
            ModuleBase::Vector3<double> displacement = velocity_for(index, step);
            for (int axis = 0; axis < 3; ++axis)
            {
                if (!reference_types[0].mbl[index][axis])
                {
                    displacement[axis] = 0.0;
                }
            }
            reference_types[0].taud[index] += displacement;
            reference_types[0].dis[index] = displacement;
        }
        unitcell::periodic_boundary_adjustment(reference_types.data(), lattice, 1);
        for (LocalAtom& atom : cell.owned_atoms())
        {
            const int index = atom.type_index;
            atom.vel = velocity_for(index, step);
        }
        integrator.first_half(log);
        for (const LocalAtom& atom : cell.owned_atoms())
        {
            expect_local_state(atom);
        }
        // Before migration, rank zero and the other nonowners can be empty.
        backing.ionic_position_updated = false;
        cell.sync_backing_unitcell(true);
        expect_backing_state(true);
        decomposition.migrate_owned_atoms(cell);
        int migrated = 0;
        for (const LocalAtom& atom : cell.owned_atoms())
        {
            expect_local_state(atom);
            const int owner = decomposition.owner_rank_from_frac(atom.frac);
            EXPECT_EQ(atom.owner_rank, owner);
            if (owner != domain.size() - 1)
            {
                ++migrated;
            }
        }
        const int global_migrated = domain.sum(migrated);
        if (domain.size() > 1)
        {
            EXPECT_GT(global_migrated, 0);
        }
        for (const LocalAtom& ghost : cell.ghost_atoms())
        {
            expect_local_state(ghost);
        }
        // Restart output may synchronize the same MD step a second time.
        for (int repeat = 0; repeat < 2; ++repeat)
        {
            const bool positions_updated = repeat == 0;
            if (positions_updated)
            {
                backing.ionic_position_updated = false;
            }
            cell.sync_backing_unitcell(positions_updated);
            expect_backing_state(true);
        }
        expect_phase_correction(false);
        expect_phase_correction(true);
        const int crossing = step == 0 ? 1 : 0;
        EXPECT_EQ(backing_types[0].boundary_shift[0].x, -crossing);
        EXPECT_EQ(backing_types[0].boundary_shift[1].x, crossing);
        EXPECT_EQ(backing_types[0].boundary_shift[2].norm2(), 0);
        EXPECT_EQ(backing_types[0].boundary_shift[3].norm2(), 0);
    }
}

TEST_F(MdBoundaryTest, HybridBoundaryMatchesFreshFoldingAtSavedFieldTime)
{
    const int rank = domain.rank();
    MD_base integrator(parameters, false, true, rank, cell);
    for (LocalAtom& atom : cell.owned_atoms())
    {
        const int index = atom.type_index;
        atom.vel = velocity_for(index, 0);
    }
    std::ofstream log;
    integrator.first_half(log);
    cell.sync_backing_unitcell(true);
    EXPECT_EQ(backing_types[0].boundary_shift[0].x, -1);
    EXPECT_EQ(backing_types[0].boundary_shift[1].x, 1);
    expect_hybrid_fold_correction(false);
    expect_hybrid_fold_correction(true);
}

TEST_F(MdBoundaryTest, InitialForceEvaluationDoesNotAdvancePositionUpdateState)
{
    const int rank = domain.rank();
    MD_base integrator(parameters, false, true, rank, cell);
    PositionStateSolver solver;
    const std::string readin_directory;
    backing.ionic_position_updated = false;
    integrator.setup(&solver, readin_directory, decomposition);
    ASSERT_EQ(solver.observed_position_updates.size(), 1u);
    EXPECT_FALSE(solver.observed_position_updates[0]);
    // The legacy setup marks the state only after the first force evaluation.
    EXPECT_TRUE(backing.ionic_position_updated);

    backing.ionic_position_updated = false;
    for (LocalAtom& atom : cell.owned_atoms())
    {
        const int index = atom.type_index;
        atom.vel = velocity_for(index, 0);
    }
    std::ofstream log;
    integrator.first_half(log);
    MD_func::force_virial(&solver, 1, cell, decomposition,
                          integrator.potential, false, integrator.virial, false);
    ASSERT_EQ(solver.observed_position_updates.size(), 2u);
    EXPECT_TRUE(solver.observed_position_updates[1]);
}

TEST_F(MdBoundaryTest, InitialSyncPreservesPositionUpdateState)
{
    for (int index = 0; index < 4; ++index)
    {
        backing_types[0].boundary_shift[index].set(7, 8, 9);
        backing_types[0].dis[index].set(7.0, 8.0, 9.0);
    }
    const bool positions_updated = false;
    for (const bool initial_state : {false, true})
    {
        backing.ionic_position_updated = initial_state;
        cell.sync_backing_unitcell(positions_updated);
        expect_backing_state(initial_state);
    }
}
} // namespace

int main(int argc, char** argv)
{
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    ::testing::InitGoogleTest(&argc, argv);
    const int local_result = RUN_ALL_TESTS();
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const int result = domain.max(local_result);
    Parallel_Global::finalize_mpi();
    return result;
}
