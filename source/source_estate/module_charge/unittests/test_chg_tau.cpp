#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_base/module_mixing/mixing.h"
#include "source_base/module_mixing/plain_mixing.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_tau.h"

#include <complex>
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
 *  unit test of module_charge/chg_tau.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - mix_tau_recip: mixes the kinetic-energy density in reciprocal space.
 *     Covered:
 *       - null pointer abort paths (chr, rhopw, rhodpw, mixing).
 *       - nspin < 1 abort.
 *       - double_grid with null mixing_highf abort.
 *       - non-double-grid plain mixing value for nspin=1.
 */

class ChgTauTest : public ::testing::Test
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

    void setup_charge(int nspin)
    {
        charge.set_rhopw(&pw_basis);
        const bool kin_den = true;
        const bool meta_gga = false;
        charge.allocate(nspin, kin_den, meta_gga, 0);
    }
};

TEST_F(ChgTauTest, NullChrAborts)
{
    Base_Mixing::Plain_Mixing mixing;
    Base_Mixing::Mixing_Data mdata;
    EXPECT_DEATH(module_charge::detail::mix_tau_recip(
                     nullptr, 1, false, &pw_basis, &pw_basis, &mixing, mdata, nullptr),
                 "");
}

TEST_F(ChgTauTest, NullGridAborts)
{
    setup_charge(1);
    Base_Mixing::Plain_Mixing mixing;
    Base_Mixing::Mixing_Data mdata;
    EXPECT_DEATH(module_charge::detail::mix_tau_recip(
                     &charge, 1, false, nullptr, &pw_basis, &mixing, mdata, nullptr),
                 "");
}

TEST_F(ChgTauTest, NullMixingAborts)
{
    setup_charge(1);
    Base_Mixing::Mixing_Data mdata;
    EXPECT_DEATH(module_charge::detail::mix_tau_recip(
                     &charge, 1, false, &pw_basis, &pw_basis, nullptr, mdata, nullptr),
                 "");
}

TEST_F(ChgTauTest, BadNspinAborts)
{
    setup_charge(1);
    Base_Mixing::Plain_Mixing mixing;
    Base_Mixing::Mixing_Data mdata;
    EXPECT_DEATH(module_charge::detail::mix_tau_recip(
                     &charge, 0, false, &pw_basis, &pw_basis, &mixing, mdata, nullptr),
                 "");
}

TEST_F(ChgTauTest, DoubleGridWithoutHighfAborts)
{
    setup_charge(1);
    Base_Mixing::Plain_Mixing mixing;
    Base_Mixing::Mixing_Data mdata;
    EXPECT_DEATH(module_charge::detail::mix_tau_recip(
                     &charge, 1, true, &pw_basis, &pw_basis, &mixing, mdata, nullptr),
                 "");
}

TEST_F(ChgTauTest, NonDoubleGridPlainMixingValue)
{
    const int nspin = 1;
    setup_charge(nspin);
    Base_Mixing::Plain_Mixing mixing(0.5);
    Base_Mixing::Mixing_Data mdata;
    mixing.init_mixing_data(mdata, pw_basis.npw, sizeof(std::complex<double>));

    // uniform kinetic densities
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.kin_r_save[0][ir] = 2.0;
        charge.kin_r[0][ir] = 3.0;
    }

    module_charge::detail::mix_tau_recip(
        &charge, nspin, false, &pw_basis, &pw_basis, &mixing, mdata, nullptr);

    // after plain mixing: out = in + beta * (out - in) = 2 + 0.5*(3-2) = 2.5
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_NEAR(charge.kin_r[0][ir], 2.5, 1e-6);
    }
}
