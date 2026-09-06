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
    Run_MD::prepare_mdcell(mdcell, ucell);

    EXPECT_EQ(mdcell.nat(), ucell.nat);
    EXPECT_EQ(mdcell.stru_meta().species.size(), 1U);
}
