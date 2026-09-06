/**********************************************
 *  Unit tests for Plus_U_Base::init_base.
 *
 *  Focus: the Yukawa state must follow the
 *  yukawa_potential argument on every call,
 *  including a true -> false re-initialization
 *  (before_scf() -> setup_pot() may call
 *  init_base() repeatedly on the same object).
 ***********************************************/

#include "source_pw/module_pwdft/dftu_base.h"

#include "source_cell/atom_spec.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

#include <vector>

class DFTUBaseTest : public testing::Test
{
  protected:
    UnitCell ucell;

    void SetUp() override
    {
        // Minimal one-atom cell: d channel available (nwl = 2),
        // one chi for each of s / p / d, so nw = 1 + 3 + 5 = 9.
        const int nw = 9;

        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.atoms[0].taud.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.itia2iat(0, iat) = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.atoms[0].taud[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        }
        ucell.atoms[0].na = 1;
        ucell.atoms[0].label = "Fe";
        ucell.atoms[0].nwl = 2;
        ucell.atoms[0].l_nchi = {1, 1, 1};
        ucell.atoms[0].nw = nw;
        ucell.atoms[0].iw2l.resize(nw);
        ucell.atoms[0].iw2n.resize(nw);
        ucell.atoms[0].iw2m.resize(nw);
        int iw = 0;
        for (int l = 0; l <= ucell.atoms[0].nwl; l++)
        {
            for (int m = 0; m < 2 * l + 1; m++)
            {
                ucell.atoms[0].iw2l[iw] = l;
                ucell.atoms[0].iw2n[iw] = 0;
                ucell.atoms[0].iw2m[iw] = m;
                iw++;
            }
        }
        ucell.set_iat2iwt(1);
    }

    void TearDown() override
    {
        // set_atom_flag is false, so ~UnitCell() skips atoms but frees
        // iat2it / iat2ia itself; only atoms must be deleted here.
        delete[] ucell.atoms;
    }

    /// Call init_base with the given Yukawa switch on a fresh d orbital
    void init_dftu(Plus_U_Base& dftu, const bool yukawa_potential)
    {
        const std::vector<int> orbital_corr = {2};
        const std::vector<double> hubbard_u = {0.0};
        dftu.init_base(ucell,
                       1,                      // npol
                       2,                      // nspin
                       orbital_corr,
                       yukawa_potential,
                       0.5,                    // yukawa_lambda
                       "",                     // global_readin_dir
                       "",                     // global_out_dir
                       "none",                 // init_chg
                       "cpu",                  // device
                       1,                      // kpar
                       hubbard_u,
                       0.0,                    // uramping
                       0,                      // occ_mat_ctrl
                       0);                     // mixing_dftu
    }
};

/// After a true -> false re-initialization the Yukawa object must be
/// released so that use_yukawa() reflects the latest argument.
TEST_F(DFTUBaseTest, InitBaseYukawaTrueThenFalseClearsState)
{
    Plus_U_Base dftu;

    init_dftu(dftu, true);
    EXPECT_TRUE(dftu.use_yukawa());

    init_dftu(dftu, false);
    EXPECT_FALSE(dftu.use_yukawa());
}

/// A false -> true re-initialization must create the Yukawa object.
TEST_F(DFTUBaseTest, InitBaseYukawaFalseThenTrueCreatesObject)
{
    Plus_U_Base dftu;

    init_dftu(dftu, false);
    EXPECT_FALSE(dftu.use_yukawa());

    init_dftu(dftu, true);
    EXPECT_TRUE(dftu.use_yukawa());
}
