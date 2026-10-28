#include "source_md/md_base.h"
#include "source_md/md_func.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_global.h"
#include "source_cell/domain_decomposition.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"
#include "source_lcao/module_rt/boundary_fix.h"

#include <gtest/gtest.h>

#include <array>
#include <complex>
#include <fstream>
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
        const double kx = kpoint == 0 ? 0.0 : 0.25;
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
        std::vector<ModuleBase::Vector3<double>> points(2);
        points[1].set(0.25, 0.0, 0.0);
        const int rows = gathered ? basis_size : orbitals.get_row_size();
        const int columns = gathered ? basis_size : orbitals.get_col_size();
        const int matrix_size = rows * columns;
        const ct::TensorShape shape({2, matrix_size});
        ct::Tensor h(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        ct::Tensor s(ct::DataType::DT_COMPLEX_DOUBLE, shape);
        const std::vector<int> plane_wave_counts;
        psi::Psi<Complex> coefficients(2, orbitals.ncol_bands, orbitals.nrow, plane_wave_counts, true);
        for (int kpoint = 0; kpoint < 2; ++kpoint)
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
            for (int element = 0; element < orbitals.nrow * orbitals.ncol_bands; ++element)
            {
                coefficients.get_pointer()[element] = Complex(1.0, 0.5);
            }
        }
        module_rt::reset_matrix_boundary(backing, points, 2, &orbitals, h, s, &coefficients, matrix_size);
        for (int kpoint = 0; kpoint < 2; ++kpoint)
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
            for (int column = 0; column < orbitals.ncol_bands; ++column)
            {
                for (int row = 0; row < orbitals.nrow; ++row)
                {
                    const int global_row = orbitals.local2global_row(row);
                    const Complex expected = Complex(1.0, 0.5) * std::conj(phase_for(global_row, kpoint));
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
