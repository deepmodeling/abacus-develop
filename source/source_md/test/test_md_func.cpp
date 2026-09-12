#include "source_md/md_func.h"
#include "source_cell/mdcell.h"
#include "source_cell/module_neighlist/domain_decomposition.h"
#include "source_base/parallel_cell.h"
#include "source_esolver/esolver.h"

#include <gtest/gtest.h>

namespace
{
class GhostForceSolver : public ModuleESolver::ESolver
{
public:
    void before_all_runners(BaseCell&, const Input_para&) override {}
    void after_all_runners(BaseCell&) override {}
    void runner(BaseCell& base, const int) override
    {
        MDCell& cell = static_cast<MDCell&>(base);
        EXPECT_TRUE(cell.has_neighbor_search());
        EXPECT_GT(cell.ghost_atoms().size(), 0);
        for (LocalAtom& atom : cell.owned_atoms())
        {
            atom.force.set(4.0, 0.0, 0.0);
        }
        for (LocalAtom& atom : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(atom.force.norm2(), 0.0);
            atom.force.set(2.0, 0.0, 0.0);
        }
    }
    double cal_energy() override { return 6.0; }
    void cal_force(BaseCell&, ModuleBase::matrix&) override
    {
        ADD_FAILURE() << "Direct MDCell forces must not use cal_force.";
    }
    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override
    {
        stress.zero_out();
        stress(0, 0) = 8.0;
    }
};

TEST(MdForceVirialTest, ReturnsGhostForcesBeforeConvertingUnitsExactlyOnce)
{
    MDCell cell;
    DomainDecomposition decomp;
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 4.0;
    lattice.e22 = 4.0;
    lattice.e33 = 4.0;
    LocalAtom atom;
    atom.frac.set(0.05, 0.5, 0.5);
    atom.cart = atom.frac * lattice;
    cell.initialize_from_owned_atoms(lattice, lattice.Inverse(), 1.0, 64.0, 1,
                                     std::vector<LocalAtom>(1, atom),
                                     std::vector<std::string>(1, "X"),
                                     std::vector<double>(1, 1.0),
                                     std::vector<std::int64_t>(1, 1), 0.4, domain);
    decomp.init(domain, lattice, 1.0, 0.0, 0.4);
    cell.set_neighbor_cutoff(0.5);
    GhostForceSolver solver;
    double potential = 0.0;
    ModuleBase::matrix virial(3, 3);
    for (int step = 0; step < 2; ++step)
    {
        MD_func::force_virial(&solver, step, cell, decomp, potential, true, virial, false);
        EXPECT_DOUBLE_EQ(potential, 3.0);
        EXPECT_DOUBLE_EQ(virial(0, 0), 4.0);
        ASSERT_EQ(cell.owned_atoms().size(), 1);
        EXPECT_DOUBLE_EQ(cell.owned_atoms()[0].force.x, 2.0 + cell.ghost_atoms().size());
        for (const LocalAtom& ghost : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(ghost.force.norm2(), 0.0);
        }
    }
}
} // namespace
