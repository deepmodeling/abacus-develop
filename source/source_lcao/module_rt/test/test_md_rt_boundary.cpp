#include "boundary_test.h"
#include "source_cell/dom_decomp.h"
#include "source_lcao/module_rt/td_folding.h"
#include "source_md/md_base.h"

#include <map>

namespace
{
using rt_test::Complex;
using rt_test::Vector;

// Fold the same HR before/after MD wrapping, keeping the original saved-A oracle.
std::vector<Complex> fold_matrix(
    const UnitCell &cell,
    const Parallel_Orbitals &layout,
    const Vector &point,
    const Vector &field,
    bool wrapped)
{
    hamilt::HContainer<Complex> hr(&layout);
    for (int i = 0; i < cell.nat; ++i)
    {
        for (int j = 0; j < cell.nat; ++j)
        {
            if (layout.get_nrow_atom(i) == 0 || layout.get_ncol_atom(j) == 0)
            {
                continue;
            }
            const int rx = j - i;
            const int ry = i - j;
            ModuleBase::Vector3<int> r(rx, ry, 0);
            if (wrapped)
            {
                r += cell.atoms[0].boundary_shift[i] - cell.atoms[0].boundary_shift[j];
            }
            hamilt::AtomPair<Complex> pair(
                i,
                j,
                r,
                &layout);
            hr.insert_pair(pair);
        }
    }
    hr.allocate(nullptr, true);
    std::map<ModuleBase::Vector3<int>, Complex> phases;
    for (int index = 0; index < hr.size_atom_pairs(); ++index)
    {
        auto &pair = hr.get_atom_pair(index);
        const auto r = pair.get_R_index(0);
        const Vector direct(r.x, r.y, r.z);
        const Vector cartesian = direct * cell.latvec * cell.lat0;
        const double angle = field * cartesian;
        const Complex exponent(0.0, angle);
        phases[r] = std::exp(exponent);
        pair.get_pointer(0)[0] = rt_test::matrix_value(pair.get_atom_i(), pair.get_atom_j());
    }
    std::vector<Complex> folded(layout.para_2d().size());
    const int rows = layout.para_2d().nrloc();
    module_rt::folding_HR_td(
        hr,
        folded.data(),
        point,
        field,
        phases,
        rows,
        1);
    return folded;
}

class MdRtBoundaryTest : public rt_test::BoundaryTest
{
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
            ASSERT_NO_FATAL_FAILURE(check_correction(no_field, rt_test::analytic_matrix));
            ASSERT_NO_FATAL_FAILURE(check_correction(saved_field, fold_matrix));
        }
    }
}

INSTANTIATE_TEST_SUITE_P(MatrixLayouts, MdRtBoundaryTest, testing::Bool());
} // namespace
