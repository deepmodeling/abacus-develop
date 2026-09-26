#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_init.h"
#include "source_io/module_restart/restart.h"

#include <complex>
#include <sstream>
#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

// chg_init.cpp references GlobalC::restart; provide a definition.
namespace GlobalC
{
Restart restart;
} // namespace GlobalC

/************************************************
 *  unit test of module_charge/chg_init.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - init_rho: the SCF charge-density initialization orchestrator.
 *     Covered:
 *       - init_chg == "wfc" with wfcpw == nullptr triggers WARNING_QUIT.
 *       - init_chg == "atomic" with ntype == 0 runs the atomic fallback
 *         (and Thomas-Fermi tau init when meta_gga is true) without crashing.
 */

namespace
{

module_charge::InitRhoCfg make_init_cfg(const std::string& init_chg, bool meta_gga)
{
    return {init_chg, "", "scf", "", 1.0, 0, 0, false, false, meta_gga, 1};
}

} // namespace

class ChgInitTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    Charge charge;
    UnitCell ucell;
    Parallel_Grid pgrid;

    void SetUp() override
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        charge.set_rhopw(&pw_basis);
        ucell.ntype = 0;
        ucell.nat = 0;
        ucell.omega = pw_basis.omega;
    }
};

TEST_F(ChgInitTest, InitChgWfcWithNullWfcpwAborts)
{
    charge.allocate(1, false, false, 0);
    ModuleSymmetry::Symmetry symm;
    ModuleBase::ComplexMatrix strucFac(0, pw_basis.npw);
    module_charge::InitRhoCfg cfg = make_init_cfg("wfc", false);

    EXPECT_DEATH(module_charge::init_rho(charge, pw_basis, ucell, pgrid, strucFac,
                                         symm, nullptr, nullptr, cfg),
                 "");
}

TEST_F(ChgInitTest, InitChgAtomicNtypeZeroMetaGgaRuns)
{
    const bool meta_gga = true;
    charge.allocate(1, meta_gga, false, 0);
    ModuleSymmetry::Symmetry symm;
    ModuleBase::ComplexMatrix strucFac(0, pw_basis.npw);
    module_charge::InitRhoCfg cfg = make_init_cfg("atomic", meta_gga);

    // ntype==0: atomic_rho loop is skipped; TF tau is computed from rho.
    module_charge::init_rho(charge, pw_basis, ucell, pgrid, strucFac,
                            symm, nullptr, nullptr, cfg);

    // kin_r should be the Thomas-Fermi expression: fact * |rho|^(5/3).
    const double fact = (3.0 / 5.0) * std::pow(3.0 * ModuleBase::PI * ModuleBase::PI, 2.0 / 3.0);
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        const double expected = fact * std::pow(std::abs(charge.rho[0][ir]), 5.0 / 3.0);
        EXPECT_NEAR(charge.kin_r[0][ir], expected, 1e-6);
    }
}

TEST_F(ChgInitTest, InitChgAutoMissingFileFallsBackToAtomic)
{
    charge.allocate(1, false, false, 0);
    ModuleSymmetry::Symmetry symm;
    ModuleBase::ComplexMatrix strucFac(0, pw_basis.npw);
    module_charge::InitRhoCfg cfg = make_init_cfg("auto", false);
    // Point the read-in directory at a path that holds no charge files.
    cfg.global_readin_dir = "no_such_dir/";

    // init_chg == "auto" must not abort when no density file exists; it
    // returns having taken the atomic-density fallback (ntype==0 keeps it a
    // no-op, but the run must survive).
    module_charge::init_rho(charge, pw_basis, ucell, pgrid, strucFac,
                            symm, nullptr, nullptr, cfg);
    SUCCEED();
}

TEST_F(ChgInitTest, InitChgFileMissingFileAborts)
{
    charge.allocate(1, false, false, 0);
    ModuleSymmetry::Symmetry symm;
    ModuleBase::ComplexMatrix strucFac(0, pw_basis.npw);
    module_charge::InitRhoCfg cfg = make_init_cfg("file", false);
    cfg.global_readin_dir = "no_such_dir/";

    // init_chg == "file" must abort when no density source can be read.
    EXPECT_DEATH(module_charge::init_rho(charge, pw_basis, ucell, pgrid, strucFac,
                                         symm, nullptr, nullptr, cfg),
                 "");
}
