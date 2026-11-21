#include "boundary_test.h"
#include "source_cell/dom_decomp.h"
#include "source_lcao/module_rt/td_folding.h"
#include "source_md/md_base.h"

namespace
{
using rt_test::Complex;
using rt_test::Vector;

class MdRtBoundaryTest : public rt_test::BoundaryTest
{
  protected:
    std::vector<Complex> fold_hybrid(
        const Parallel_Orbitals &orbitals,
        const Vector &kpoint,
        const Vector &vector_potential,
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
                    translation += types[0].boundary_shift[row] - types[0].boundary_shift[column];
                }
                hamilt::AtomPair<Complex> pair(
                    row,
                    column,
                    translation,
                    &orbitals);
                real_space.insert_pair(pair);
            }
        }
        real_space.allocate(nullptr, true);
        std::map<ModuleBase::Vector3<int>, Complex> phases;
        for (int index = 0; index < real_space.size_atom_pairs(); ++index)
        {
            hamilt::AtomPair<Complex> &pair = real_space.get_atom_pair(index);
            const ModuleBase::Vector3<int> translation = pair.get_R_index(0);
            const int row = pair.get_atom_i();
            const int column = pair.get_atom_j();
            pair.get_pointer(0)[0] = rt_test::matrix_value(row, column);
            const Vector direct(translation.x, translation.y, translation.z);
            const Vector cartesian = direct * cell.latvec * cell.lat0;
            const double angle = vector_potential * cartesian;
            const Complex exponent(0.0, angle);
            phases[translation] = std::exp(exponent);
        }
        std::vector<Complex> folded(orbitals.para_2d().size(), Complex(0.0, 0.0));
        const int rows = orbitals.para_2d().nrloc();
        module_rt::folding_HR_td(
            real_space,
            folded.data(),
            kpoint,
            vector_potential,
            phases,
            rows,
            1);
        return folded;
    }

    void expect_hybrid_fold_correction(const Vector &saved_A)
    {
        const auto &layout = matrix_layout();
        const int matrix_size = layout.para_2d().size();
        const Vector next_A(-0.13, 0.17, 0.09);
        for (int k = 0; k < 4; ++k)
        {
            const auto before = fold_hybrid(
                layout,
                points[k],
                saved_A,
                false);
            const auto after = fold_hybrid(
                layout,
                points[k],
                saved_A,
                true);
            const auto wrong_time = fold_hybrid(
                layout,
                points[k],
                next_A,
                true);
            double field_time_difference = 0.0;
            for (int element = 0; element < matrix_size; ++element)
            {
                const int offset = k * matrix_size + element;
                h.data<Complex>()[offset] = before[element];
                expected[offset] = after[element];
                field_time_difference += std::abs(after[element] - wrong_time[element]);
            }
            EXPECT_GT(field_time_difference, 0.01);
        }
        expect_correction(saved_A);
    }
};

TEST_P(MdRtBoundaryTest, SavedStateSurvivesWrappingMigrationAndRepeatedSync)
{
    const std::array<Vector, 4> positions = {
        {{0.99, 0.25, 0.75}, {0.01, 0.75, 0.25}, {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}};
    Atom &type = types[0];
    type.tau.resize(4);
    type.taud.resize(4);
    type.dis.resize(4);
    type.vel.resize(4);
    type.mbl.resize(4, ModuleBase::Vector3<int>(1, 1, 1));
    type.mbl[3].set(0, 0, 0);
    std::vector<LocalAtom> owned;
    const int source_rank = domain.size() - 1;
    for (int i = 0; i < 4; ++i)
    {
        type.taud[i] = positions[i];
        type.tau[i] = positions[i] * cell.latvec;
        if (domain.rank() == source_rank)
        {
            LocalAtom atom;
            atom.frac = type.taud[i];
            atom.cart = type.tau[i];
            atom.mbl = type.mbl[i];
            atom.type_index = i;
            atom.owner_rank = source_rank;
            owned.push_back(atom);
        }
    }
    MDCell mdcell;
    const auto inverse = cell.latvec.Inverse();
    const std::vector<std::string> labels(1, "X");
    const std::vector<double> masses(1, 1.0);
    const std::vector<std::int64_t> counts(1, 4);
    const double omega = cell.lat0 * cell.lat0 * cell.lat0;
    mdcell.initialize_from_owned_atoms(
        cell.latvec,
        inverse,
        cell.lat0,
        omega,
        4,
        owned,
        labels,
        masses,
        counts,
        0.0,
        domain);
    mdcell.set_backing_unitcell(cell);
    mdcell.set_neighbor_cutoff(0.05);
    DomainDecomposition decomposition;
    decomposition.init(
        domain,
        cell.latvec,
        cell.lat0,
        0.0,
        0.0);
    MD_para parameters;
    parameters.md_dt = ModuleBase::AU_to_FS;
    parameters.md_tfirst = -1.0;
    parameters.md_tlast = -1.0;
    const int rank = domain.rank();
    MD_base integrator(
        parameters,
        false,
        true,
        rank,
        mdcell);
    std::ofstream log;
    const Vector no_field;
    const Vector saved_field(0.07, -0.11, 0.05);
    for (int step = 0; step < 2; ++step)
    {
        SCOPED_TRACE(step);
        for (LocalAtom &atom : mdcell.owned_atoms())
        {
            Vector drift(0.0, 0.01, 0.0);
            if (step == 0)
            {
                drift.set(-1.0e-21, -2.0e-21, -3.0e-21);
                if (atom.type_index == 0)
                {
                    drift.set(0.02, 0.6, 0.0);
                }
                if (atom.type_index == 1)
                {
                    drift.set(-0.02, -0.6, 0.0);
                }
            }
            atom.vel = drift * cell.latvec * cell.lat0;
        }
        integrator.first_half(log);
        mdcell.sync_backing_unitcell(true);
        if (step == 0)
        {
            ASSERT_NO_FATAL_FAILURE(expect_hybrid_fold_correction(saved_field));
        }
        decomposition.migrate_owned_atoms(mdcell);
        for (int repeat = 0; repeat < 2; ++repeat)
        {
            const bool positions_updated = repeat == 0;
            mdcell.sync_backing_unitcell(positions_updated);
            const int crossing = step == 0 ? 1 : 0;
            EXPECT_EQ(type.boundary_shift[0].x, -crossing);
            EXPECT_EQ(type.boundary_shift[1].x, crossing);
            EXPECT_EQ(type.boundary_shift[2].norm2(), 0);
            EXPECT_EQ(type.boundary_shift[3].norm2(), 0);
        }
        ASSERT_NO_FATAL_FAILURE(expect_phase_correction(no_field));
    }
}

INSTANTIATE_TEST_SUITE_P(MatrixLayouts, MdRtBoundaryTest, testing::Bool());
TEST(TdFoldingTest, BlochAndHybridPhases)
{
    using Complex = std::complex<double>;
    using Vector = ModuleBase::Vector3<double>;
    Parallel_Orbitals orbitals;
    orbitals.para_2d().set_serial(1, 1);
    const int atom_start = 0;
    orbitals.set_atomic_trace(&atom_start, 1, 1);
    const ModuleBase::Vector3<int> r(1, -1, 2);
    hamilt::HContainer<Complex> hr(&orbitals);
    hamilt::AtomPair<Complex> pair(
        0,
        0,
        r,
        &orbitals);
    hr.insert_pair(pair);
    hr.allocate(nullptr, true);
    const Complex value(1.0, 0.2);
    hr.get_atom_pair(0).get_pointer(0)[0] = value;
    const Vector k(0.25, -0.125, 0.375);
    const Vector direct(r.x, r.y, r.z);
    const std::vector<Vector> fields = {{0.0, 0.0, 0.0}, {0.07, -0.11, 0.05}};
    for (const auto &field : fields)
    {
        const double field_angle = field * direct;
        const Complex field_exponent(0.0, field_angle);
        const std::map<ModuleBase::Vector3<int>, Complex> phases = {{r, std::exp(field_exponent)}};
        Complex folded(0.0, 0.0);
        module_rt::folding_HR_td(
            hr,
            &folded,
            k,
            field,
            phases,
            1,
            1);
        const double angle = ModuleBase::TWO_PI * (k * direct) + field_angle;
        const Complex exponent(0.0, angle);
        const Complex expected = value * std::exp(exponent);
        const double error = std::abs(folded - expected);
        EXPECT_LT(error, 1.0e-12);
    }
}

} // namespace
