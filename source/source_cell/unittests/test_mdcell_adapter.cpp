#include "source_cell/mdcell_adapter.h"

#include "source_base/matrix.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_global.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

TEST(MDCellAdapterTest, SyncGeometryAndCollectForces)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.lat0 = 1.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e33 = 1.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.omega = 1.0;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Ar";
    ucell.atoms[0].mass = 39.948;
    ucell.atoms[0].na = 2;
    ucell.atoms[0].tau.resize(2);
    ucell.atoms[0].taud.resize(2);
    ucell.atoms[0].vel.resize(2);
    ucell.atoms[0].mbl.resize(2);
    ucell.atoms[0].tau[0].set(0.1, 0.1, 0.1);
    ucell.atoms[0].tau[1].set(0.7, 0.7, 0.7);

    MDCell mdcell;
    DomainDecomposition decomp;
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    MDCellAdapter::initialize_from_unitcell(mdcell, decomp, ucell, 0.0, domain);
    mdcell.set_neighbor_cutoff(0.4);

    ucell.lat0 = 1.5;
    ucell.latvec.e11 = 2.0;
    ucell.latvec.e22 = 2.0;
    ucell.latvec.e33 = 2.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.omega = 27.0;
    ucell.atoms[0].tau[0].set(0.3, 0.4, 0.5);
    ucell.atoms[0].tau[1].set(1.4, 1.5, 1.6);
    MDCellAdapter::sync_geometry_from_unitcell(mdcell, ucell);
    EXPECT_DOUBLE_EQ(mdcell.lat0(), 1.5);
    EXPECT_DOUBLE_EQ(mdcell.omega(), 27.0);
    for (const LocalAtom& atom : mdcell.owned_atoms())
    {
        const double expected_x = atom.type_index == 0 ? 0.3 : 1.4;
        EXPECT_NEAR(atom.cart.x, expected_x, 1.0e-12);
    }

    decomp.prepare_neighbors(mdcell);
    for (LocalAtom& atom : mdcell.owned_atoms())
    {
        const double force = static_cast<double>(atom.type_index + 1);
        atom.force.set(force, 2.0 * force, 3.0 * force);
    }
    ModuleBase::matrix full_force;
    MDCellAdapter::collect_full_forces(mdcell, decomp, ucell, full_force);
    EXPECT_DOUBLE_EQ(full_force(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(full_force(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(full_force(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(full_force(1, 2), 6.0);
}

int main(int argc, char** argv)
{
    int nproc = 1;
    int nthreads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, nproc, nthreads, rank);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
