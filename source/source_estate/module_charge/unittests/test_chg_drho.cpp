#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_drho.h"
#include "source_estate/module_charge/chg_mix_cfg.h"

#include <cmath>
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

/************************************************
 *  unit test of module_charge/chg_drho.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - inner_product_real: sum of rho1*rho2 over nrxx*nspin/resize_tmp
 *     - nspin == 1: single block
 *     - nspin == 4 with mixing_angle > 0: resize_tmp == 2 (half length)
 *   - cal_drho: real-space (scf_thr_type == 2) L1 residual normalized by nelec
 *     - nspin == 1: sum over single channel
 *     - nspin == 2: both channels
 *     - nspin == 4 with domag_z: only channels 0 and 3
 *   - cal_dkin: kinetic-energy-density residual
 *     - meta_gga == false: returns 0
 *     - meta_gga == true: same L1 logic as cal_drho real-space
 */

namespace
{

MixingConfig make_cfg(int nspin, int scf_thr_type, bool domag_z = false, double mixing_angle = -10.0)
{
    MixingConfig cfg{
        "plain",       // mixing_mode
        0.8,           // mixing_beta
        4,             // mixing_ndim
        0.0,           // mixing_gg0
        false,         // mixing_tau
        1.6,           // mixing_beta_mag
        0.0,           // mixing_gg0_mag
        0.1,           // mixing_gg0_min
        mixing_angle,  // mixing_angle
        false,         // mixing_dmr
        nspin,         // nspin
        scf_thr_type,  // scf_thr_type
        false,         // double_grid
        false,         // gamma_only_pw
        false,         // domag
        domag_z,       // domag_z
        100            // scf_nmax
    };
    return cfg;
}

} // namespace

class ChgDrhoTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    Charge charge;

    void SetUp() override
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
    }

    /// Allocate Charge buffers for nspin with kin_r allocated.
    void setup_charge(int nspin)
    {
        charge.set_rhopw(&pw_basis);
        const bool kin_den = true;
        const bool meta_gga = true;
        charge.allocate(nspin, kin_den, meta_gga, 0);
    }
};

// ---------------------------------------------------------------------------
// inner_product_real
// ---------------------------------------------------------------------------

TEST_F(ChgDrhoTest, InnerProductRealNspin1)
{
    const int nspin = 1;
    MixingConfig cfg = make_cfg(nspin, 2);
    std::vector<double> rho1(pw_basis.nrxx, 2.0);
    std::vector<double> rho2(pw_basis.nrxx, 3.0);

    const double inner = module_charge::inner_product_real(
        rho1.data(), rho2.data(), pw_basis, cfg);

    // sum of 2.0 * 3.0 over nrxx elements
    EXPECT_NEAR(inner, 6.0 * pw_basis.nrxx, 1e-8);
}

TEST_F(ChgDrhoTest, InnerProductRealNspin4AngleHalvesLength)
{
    const int nspin = 4;
    MixingConfig cfg = make_cfg(nspin, 2, false, 1.0); // mixing_angle > 0
    const int len = pw_basis.nrxx * nspin / 2; // resize_tmp == 2
    std::vector<double> rho1(len, 1.0);
    std::vector<double> rho2(len, 1.0);

    const double inner = module_charge::inner_product_real(
        rho1.data(), rho2.data(), pw_basis, cfg);

    EXPECT_NEAR(inner, static_cast<double>(len), 1e-8);
}

// ---------------------------------------------------------------------------
// cal_drho (real-space path, scf_thr_type == 2)
// ---------------------------------------------------------------------------

TEST_F(ChgDrhoTest, CalDrhoRealNspin1)
{
    const int nspin = 1;
    setup_charge(nspin);
    MixingConfig cfg = make_cfg(nspin, 2);

    const double nelec = 4.0;
    const double omega = 1.0;

    // rho - rho_save = 1.0 everywhere
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.rho[0][ir] = 2.0;
        charge.rho_save[0][ir] = 1.0;
    }

    const double drho = module_charge::cal_drho(
        &charge, nelec, pw_basis, cfg, omega, 1.0);

    // drho = sum|diff| * omega/nxyz / nelec = nrxx * 1.0 * 1.0/nxyz / 4.0
    const double ref = static_cast<double>(pw_basis.nrxx) * 1.0
                       * omega / static_cast<double>(pw_basis.nxyz) / nelec;
    EXPECT_NEAR(drho, ref, 1e-10);
}

TEST_F(ChgDrhoTest, CalDrhoRealNspin2BothChannels)
{
    const int nspin = 2;
    setup_charge(nspin);
    MixingConfig cfg = make_cfg(nspin, 2, false);

    const double nelec = 4.0;
    const double omega = 1.0;

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_basis.nrxx; ++ir)
        {
            charge.rho[is][ir] = 3.0;
            charge.rho_save[is][ir] = 1.0; // diff = 2.0
        }
    }

    const double drho = module_charge::cal_drho(
        &charge, nelec, pw_basis, cfg, omega, 1.0);

    // both channels contribute: 2 * nrxx * 2.0 * omega/nxyz / nelec
    const double ref = 2.0 * pw_basis.nrxx * 2.0
                       * omega / static_cast<double>(pw_basis.nxyz) / nelec;
    EXPECT_NEAR(drho, ref, 1e-10);
}

TEST_F(ChgDrhoTest, CalDrhoRealNspin4DomagZOnlyChannels0And3)
{
    const int nspin = 4;
    setup_charge(nspin);
    MixingConfig cfg = make_cfg(nspin, 2, true); // domag_z == true

    const double nelec = 4.0;
    const double omega = 1.0;

    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw_basis.nrxx; ++ir)
        {
            charge.rho[is][ir] = 3.0;
            charge.rho_save[is][ir] = 1.0; // diff = 2.0
        }
    }

    const double drho = module_charge::cal_drho(
        &charge, nelec, pw_basis, cfg, omega, 1.0);

    // only is==0 and is==3 are summed when domag_z is true
    const double ref = 2.0 * pw_basis.nrxx * 2.0
                       * omega / static_cast<double>(pw_basis.nxyz) / nelec;
    EXPECT_NEAR(drho, ref, 1e-10);
}

// ---------------------------------------------------------------------------
// cal_dkin
// ---------------------------------------------------------------------------

TEST_F(ChgDrhoTest, CalDkinMetaGgaFalseReturnsZero)
{
    const int nspin = 1;
    setup_charge(nspin);
    charge.meta_gga = false;
    MixingConfig cfg = make_cfg(nspin, 2);

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.kin_r[0][ir] = 5.0;
        charge.kin_r_save[0][ir] = 1.0;
    }

    const double dkin = module_charge::cal_dkin(
        &charge, 4.0, pw_basis, cfg, 1.0);

    EXPECT_NEAR(dkin, 0.0, 1e-12);
}

TEST_F(ChgDrhoTest, CalDkinMetaGgaTrueComputesResidual)
{
    const int nspin = 1;
    setup_charge(nspin);
    charge.meta_gga = true;
    MixingConfig cfg = make_cfg(nspin, 2);

    const double nelec = 4.0;
    const double omega = 1.0;

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.kin_r[0][ir] = 4.0;
        charge.kin_r_save[0][ir] = 1.0; // diff = 3.0
    }

    const double dkin = module_charge::cal_dkin(
        &charge, nelec, pw_basis, cfg, omega);

    const double ref = static_cast<double>(pw_basis.nrxx) * 3.0
                       * omega / static_cast<double>(pw_basis.nxyz) / nelec;
    EXPECT_NEAR(dkin, ref, 1e-10);
}
