#include "source_cell/domain_decomposition.h"
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

    EXPECT_EQ(mdcell.nat(), ucell.nat);
    EXPECT_EQ(mdcell.stru_meta().species.size(), 1U);

    // UnitCell-backed initialization uses the same wrapping as distributed input.
    ucell.atoms[0].tau[0].set(-0.2, 1.25, 2.5);
    ucell.atoms[0].taud[0] = ucell.atoms[0].tau[0];
    ucell.atoms[0].vel[0].set(1.0, 2.0, 3.0);
    Run_MD::prepare_mdcell(mdcell, ucell, decomp);
    ASSERT_EQ(mdcell.owned_atoms().size(), 1);
    EXPECT_TRUE(mdcell.has_backing_unitcell());
    EXPECT_EQ(&mdcell.backing_unitcell(), &ucell);
    const LocalAtom& atom = mdcell.owned_atoms()[0];
    EXPECT_NEAR(atom.frac.x, 0.8, 1.0e-12);
    EXPECT_NEAR(atom.frac.y, 0.25, 1.0e-12);
    EXPECT_NEAR(atom.frac.z, 0.5, 1.0e-12);
    EXPECT_DOUBLE_EQ((atom.cart - atom.frac * mdcell.latvec()).norm2(), 0.0);
    EXPECT_DOUBLE_EQ(atom.vel.y, 2.0);
    EXPECT_EQ(atom.mbl.x, 0);
    EXPECT_EQ(atom.owner_rank, 0);
}

// --- Cadence independence: output (out_freq_ion) and restart
// (md_restartfreq) must be driven by their own frequencies so that
// setting out_freq_ion does not silently re-aim md_restartfreq.

TEST(RunMDCadenceTest, disabled_when_freq_non_positive)
{
    for (int step = 0; step <= 20; ++step)
    {
        EXPECT_FALSE(Run_MD::is_cadence_step(step, 0));
        EXPECT_FALSE(Run_MD::is_cadence_step(step, -1));
    }
}

TEST(RunMDCadenceTest, fires_every_freq_steps)
{
    const int freq = 3;
    for (int step = 0; step <= 20; ++step)
    {
        const bool expected = (step % freq == 0);
        EXPECT_EQ(Run_MD::is_cadence_step(step, freq), expected) << "step=" << step;
    }
}

// Core regression guard: with out_freq_ion and md_restartfreq set to
// different positive values, the two cadences must disagree on at least
// one step. This guarantees md_restartfreq is no longer swallowed by
// out_freq_ion > 0.
TEST(RunMDCadenceTest, output_and_restart_are_independent)
{
    const int out_freq_ion = 4;
    const int md_restartfreq = 6;

    bool found_disagreement = false;
    for (int step = 0; step <= 60; ++step)
    {
        const bool output_step = Run_MD::is_cadence_step(step, out_freq_ion);
        const bool restart_step = Run_MD::is_cadence_step(step, md_restartfreq);
        if (output_step != restart_step)
        {
            found_disagreement = true;
            break;
        }
    }
    EXPECT_TRUE(found_disagreement)
        << "out_freq_ion=" << out_freq_ion << " and md_restartfreq=" << md_restartfreq
        << " produced identical cadences; md_restartfreq is still being ignored.";
}
