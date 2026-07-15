#include <gtest/gtest.h>

#include "source_md/md_state_view.h"
#include "source_cell/md_cell.h"

#include <mpi.h>
#include <cmath>
#include <vector>

namespace
{
void ensure_mpi_initialized()
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        int provided = 0;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
    }
}

double cell_volume(const ModuleBase::Matrix3& latvec)
{
    const double cx = latvec.e22 * latvec.e33 - latvec.e23 * latvec.e32;
    const double cy = latvec.e23 * latvec.e31 - latvec.e21 * latvec.e33;
    const double cz = latvec.e21 * latvec.e32 - latvec.e22 * latvec.e31;
    return std::abs(latvec.e11 * cx + latvec.e12 * cy + latvec.e13 * cz);
}

ModuleBase::Matrix3 diagonal_lattice()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = 0.0;
    latvec.e22 = 10.0;
    latvec.e23 = 0.0;
    latvec.e31 = 0.0;
    latvec.e32 = 0.0;
    latvec.e33 = 10.0;
    return latvec;
}
}

TEST(MdStateViewMdCellTest, MdCellViewMapsOwnedAtomStorageWithoutCopy)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = diagonal_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{
        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
        ModuleBase::Vector3<double>(5.2, 5.2, 0.0),
        ModuleBase::Vector3<double>(5.1, 0.0, 5.0),
        ModuleBase::Vector3<double>(0.0, 5.3, 5.0)};
    MdCell mdcell(latvec, latvec.Inverse(), 1.0, cell_volume(latvec), positions, MPI_COMM_SELF, 8.5, 0.0);

    ASSERT_EQ(mdcell.nlocal(), 4);
    mdcell.mutable_owned_atoms()[0].mass = 39.948;
    mdcell.mutable_owned_atoms()[1].mass = 39.948;
    mdcell.mutable_owned_atoms()[2].mass = 39.948;
    mdcell.mutable_owned_atoms()[3].mass = 39.948;
    mdcell.mutable_owned_atoms()[0].vel.set(-0.1, 0.1, 0.0);
    mdcell.mutable_owned_atoms()[1].vel.set(0.2, -0.2, 0.1);
    mdcell.mutable_owned_atoms()[2].vel.set(-0.3, 0.0, -0.1);
    mdcell.mutable_owned_atoms()[3].vel.set(0.4, 0.2, -0.2);
    mdcell.mutable_owned_atoms()[1].force.set(1.0, 2.0, 3.0);

    MdStateView view = MdStateView::from_mdcell(mdcell);

    ASSERT_EQ(view.size(), 4);
    EXPECT_TRUE(view.uses_distributed_storage());
    EXPECT_EQ(view.global_id(1), 1);
    EXPECT_EQ(view.type(1), 0);
    EXPECT_EQ(view.type_index(1), 1);
    EXPECT_DOUBLE_EQ(view.mass(1), 39.948);
    EXPECT_DOUBLE_EQ(view.cart(2).x, mdcell.owned_atoms()[2].cart.x);
    EXPECT_DOUBLE_EQ(view.frac(2).z, mdcell.owned_atoms()[2].frac.z);
    EXPECT_DOUBLE_EQ(view.vel(3).y, mdcell.owned_atoms()[3].vel.y);
    EXPECT_EQ(view.mbl(0).z, 1);
    EXPECT_DOUBLE_EQ(view.force(1).y, 2.0);

    view.cart(0).x = 0.125;
    view.vel(1).z = -6.0;
    view.mbl(2).x = 0;
    view.force(3).z = 8.0;

    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].cart.x, 0.125);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[1].vel.z, -6.0);
    EXPECT_EQ(mdcell.owned_atoms()[2].mbl.x, 0);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[3].force.z, 8.0);
}
