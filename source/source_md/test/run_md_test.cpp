#include "source_cell/module_neighlist/domain_decomposition.h"
#include "gtest/gtest.h"

#include "source_cell/mdcell.h"
#include "source_cell/unitcell.h"
#include "source_md/run_md.h"

TEST(RunMDTest, prepare_mdcell_from_unitcell)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Ar";
    ucell.atoms[0].mass = 39.948;
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau.resize(1);
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].vel.resize(1);
    ucell.atoms[0].mbl.resize(1);
    ucell.atoms[0].tau[0].set(0.0, 0.0, 0.0);
    ucell.atoms[0].taud[0].set(0.0, 0.0, 0.0);
    ucell.atoms[0].vel[0].set(0.0, 0.0, 0.0);
    ucell.atoms[0].mbl[0].set(0, 0, 0);

    MDCell mdcell;
    DomainDecomposition decomp;
    Run_MD::prepare_mdcell(mdcell, ucell, decomp);

    EXPECT_EQ(mdcell.nat_, ucell.nat);
    EXPECT_EQ(mdcell.stru_meta_.species.size(), 1U);

    // UnitCell-backed initialization uses the same wrapping as distributed input.
    ucell.atoms[0].tau[0].set(-0.2, 1.25, 2.5);
    ucell.atoms[0].taud[0] = ucell.atoms[0].tau[0];
    ucell.atoms[0].vel[0].set(1.0, 2.0, 3.0);
    Run_MD::prepare_mdcell(mdcell, ucell, decomp);
    ASSERT_EQ(mdcell.owned_atoms_.size(), 1);
    EXPECT_TRUE(mdcell.backing_unitcell_ != nullptr);
    EXPECT_EQ(mdcell.backing_unitcell_, &ucell);
    const LocalAtom& atom = mdcell.owned_atoms_[0];
    EXPECT_NEAR(atom.frac.x, 0.8, 1.0e-12);
    EXPECT_NEAR(atom.frac.y, 0.25, 1.0e-12);
    EXPECT_NEAR(atom.frac.z, 0.5, 1.0e-12);
    EXPECT_DOUBLE_EQ((atom.cart - atom.frac * mdcell.latvec_).norm2(), 0.0);
    EXPECT_DOUBLE_EQ(atom.vel.y, 2.0);
    EXPECT_EQ(atom.mbl.x, 0);
    EXPECT_EQ(atom.owner_rank, 0);
}
