#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_pot/pot_rvv10.h"
#include "source_estate/module_pot/pot_xc.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

namespace
{
class Rvv10Potential : public testing::Test
{
  protected:
    ModulePW::PW_Basis pw{"cpu", "double"};
    UnitCell cell;
    Charge charge;

    void SetUp() override
    {
        XC_Functional::set_xc_type("GGA_X_RPW86+GGA_C_PBE");
        pw.initgrids(10.0, ModuleBase::Matrix3(1, 0, 0, 0, 1.25, 0, 0, 0, 1.5), 8, 10, 12);
        pw.initparameters(false, 4.0, 1, false);
        pw.setuptransform();
        pw.collect_local_pw();
        cell.omega = pw.omega;
        cell.tpiba = pw.tpiba;
        charge.set_rhopw(&pw);
        charge.allocate(1, false, false, 0);
        const double pi = std::acos(-1.0);
        for (int x = 0; x < pw.nx; ++x)
            for (int y = 0; y < pw.ny; ++y)
                for (int z = 0; z < pw.nz; ++z)
                {
                    const int i = (x * pw.ny + y) * pw.nz + z;
                    charge.rho[0][i]
                        = 0.008 + 0.001 * std::cos(2 * pi * x / pw.nx) + 0.0007 * std::sin(2 * pi * y / pw.ny);
                    charge.rho_core[i] = 0.003 + 0.0002 * std::cos(2 * pi * z / pw.nz);
                }
    }

    void TearDown() override
    {
        XC_Functional::set_xc_type("PBE");
    }
};

// Catches replacement instead of addition, accumulated energies, omitted core in
// the energy, and contraction of vtxc with total instead of valence density.
TEST_F(Rvv10Potential, FullEnergyDerivativeAndValenceContractionWithFixedCore)
{
    double energy = 0;
    double vtxc = 0;
    ModuleBase::matrix v(1, pw.nrxx);
    ModuleBase::matrix vofk;
    elecstate::PotXC semilocal(&pw, &energy, &vtxc, &vofk);
    semilocal.cal_v_eff(&charge, &cell, v);
    elecstate::PotRvv10 potential(&pw, &energy, &vtxc);
    potential.cal_v_eff(&charge, &cell, v);
    const double original_energy = energy;
    const double original_vtxc = vtxc;
    const double dv = pw.omega / pw.nxyz;
    double contracted = 0;
    double core_contracted = 0;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        contracted += dv * charge.rho[0][i] * v(0, i);
        core_contracted += dv * charge.rho_core[i] * v(0, i);
    }
    EXPECT_NEAR(vtxc, contracted, 1.e-12);
    ASSERT_GT(std::abs(core_contracted), 1.e-3);

    ModuleBase::matrix accumulated(1, pw.nrxx);
    for (int i = 0; i < pw.nrxx; ++i)
        accumulated(0, i) = 0.37;
    semilocal.cal_v_eff(&charge, &cell, accumulated);
    potential.cal_v_eff(&charge, &cell, accumulated);
    EXPECT_DOUBLE_EQ(energy, original_energy);
    EXPECT_DOUBLE_EQ(vtxc, original_vtxc);
    for (int i = 0; i < pw.nrxx; ++i)
        EXPECT_NEAR(accumulated(0, i), 0.37 + v(0, i), 1.e-13);

    std::vector<double> initial(charge.rho[0], charge.rho[0] + pw.nrxx);
    std::vector<double> direction(pw.nrxx);
    double derivative = 0;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        direction[i] = initial[i] - 0.007;
        derivative += dv * direction[i] * v(0, i);
    }
    for (double h: {1.e-3, 3.e-4, 1.e-4})
    {
        for (int i = 0; i < pw.nrxx; ++i)
            charge.rho[0][i] = initial[i] + h * direction[i];
        ModuleBase::matrix plus(1, pw.nrxx);
        semilocal.cal_v_eff(&charge, &cell, plus);
        potential.cal_v_eff(&charge, &cell, plus);
        const double ep = energy;
        for (int i = 0; i < pw.nrxx; ++i)
            charge.rho[0][i] = initial[i] - h * direction[i];
        ModuleBase::matrix minus(1, pw.nrxx);
        semilocal.cal_v_eff(&charge, &cell, minus);
        potential.cal_v_eff(&charge, &cell, minus);
        const double em = energy;
        EXPECT_NEAR((ep - em) / (2 * h), derivative, 2.e-7) << "step=" << h;
    }
}

TEST_F(Rvv10Potential, NonlocalComponentOnlyAddsToExistingResults)
{
    std::vector<double> valence(charge.rho[0], charge.rho[0] + pw.nrxx);
    std::vector<double> total = valence;
    for (int i = 0; i < pw.nrxx; ++i)
        total[i] += charge.rho_core[i];

    const Rvv10::Evaluation nonlocal
        = Rvv10::SerialEvaluator(Rvv10::rvv10_b_default, Rvv10::rvv10_c_default).evaluate(pw, total, valence);
    double energy = 1.25;
    double vtxc = -0.75;
    ModuleBase::matrix potential(1, pw.nrxx);
    for (int i = 0; i < pw.nrxx; ++i)
        potential(0, i) = 0.37;

    elecstate::PotRvv10 component(&pw, &energy, &vtxc);
    component.cal_v_eff(&charge, &cell, potential);

    EXPECT_NEAR(energy, 1.25 + nonlocal.energy, 1.e-12);
    EXPECT_NEAR(vtxc, -0.75 + nonlocal.vtxc, 1.e-12);
    for (int i = 0; i < pw.nrxx; ++i)
        EXPECT_NEAR(potential(0, i), 0.37 + nonlocal.potential[i], 1.e-12);
}

TEST_F(Rvv10Potential, SpinPolarizedNonlocalPotentialIsSharedByBothChannels)
{
    Charge spin_charge;
    spin_charge.set_rhopw(&pw);
    spin_charge.allocate(2, false, false, 0);
    for (int i = 0; i < pw.nrxx; ++i)
    {
        spin_charge.rho[0][i] = 0.005 + 0.0006 * std::cos(0.1 * i);
        spin_charge.rho[1][i] = 0.003 + 0.0004 * std::sin(0.1 * i);
        spin_charge.rho_core[i] = charge.rho_core[i];
    }

    std::vector<double> valence(pw.nrxx, 0.0);
    for (int i = 0; i < pw.nrxx; ++i)
        valence[i] = spin_charge.rho[0][i] + spin_charge.rho[1][i];
    std::vector<double> total = valence;
    for (int i = 0; i < pw.nrxx; ++i)
        total[i] += spin_charge.rho_core[i];

    const Rvv10::Evaluation nonlocal
        = Rvv10::SerialEvaluator(Rvv10::rvv10_b_default, Rvv10::rvv10_c_default).evaluate(pw, total, valence);
    double energy = 0;
    double vtxc = 0;
    ModuleBase::matrix potential(2, pw.nrxx);
    for (int is = 0; is < 2; ++is)
        for (int i = 0; i < pw.nrxx; ++i)
            potential(is, i) = 0.0;

    elecstate::PotRvv10 component(&pw, &energy, &vtxc);
    component.cal_v_eff(&spin_charge, &cell, potential);

    EXPECT_NEAR(energy, nonlocal.energy, 1.e-12);
    EXPECT_NEAR(vtxc, nonlocal.vtxc, 1.e-12);
    for (int i = 0; i < pw.nrxx; ++i)
    {
        EXPECT_NEAR(potential(0, i), nonlocal.potential[i], 1.e-12);
        EXPECT_NEAR(potential(1, i), nonlocal.potential[i], 1.e-12);
    }
}

// The production adapters should be exactly the canonical semilocal PotXC path
// followed by the independent additive nonlocal evaluator.
TEST_F(Rvv10Potential, MatchesCanonicalSemilocalPlusNonlocalAdapters)
{
    std::vector<double> valence(charge.rho[0], charge.rho[0] + pw.nrxx);
    std::vector<double> total = valence;
    for (int i = 0; i < pw.nrxx; ++i)
        total[i] += charge.rho_core[i];

    double semilocal_energy = 0;
    double semilocal_vtxc = 0;
    ModuleBase::matrix vofk;
    elecstate::PotXC semilocal(&pw, &semilocal_energy, &semilocal_vtxc, &vofk);
    ModuleBase::matrix semilocal_v(1, pw.nrxx);
    semilocal.cal_v_eff(&charge, &cell, semilocal_v);

    double rvv10_energy = semilocal_energy;
    double rvv10_vtxc = semilocal_vtxc;
    ModuleBase::matrix rvv10_v = semilocal_v;
    elecstate::PotRvv10 rvv10(&pw, &rvv10_energy, &rvv10_vtxc);
    rvv10.cal_v_eff(&charge, &cell, rvv10_v);

    const Rvv10::Evaluation nonlocal
        = Rvv10::SerialEvaluator(Rvv10::rvv10_b_default, Rvv10::rvv10_c_default).evaluate(pw, total, valence);
    EXPECT_NEAR(rvv10_energy, semilocal_energy + nonlocal.energy, 1.e-12);
    EXPECT_NEAR(rvv10_vtxc, semilocal_vtxc + nonlocal.vtxc, 1.e-12);
    for (int i = 0; i < pw.nrxx; ++i)
        EXPECT_NEAR(rvv10_v(0, i), semilocal_v(0, i) + nonlocal.potential[i], 1.e-12);
}

// Reuses the actual adapters in one process. Selecting a different semilocal
// functional must not change the independent nonlocal component.
TEST_F(Rvv10Potential, SwitchingSemilocalFunctionalKeepsComponentsIndependent)
{
    double pbe_energy = 0;
    double pbe_vtxc = 0;
    ModuleBase::matrix pbe_vofk;
    elecstate::PotXC pbe(&pw, &pbe_energy, &pbe_vtxc, &pbe_vofk);
    XC_Functional::set_xc_type("PBE");
    ModuleBase::matrix pbe_potential(1, pw.nrxx);
    pbe.cal_v_eff(&charge, &cell, pbe_potential);
    const double reference_pbe_energy = pbe_energy;
    const double reference_pbe_vtxc = pbe_vtxc;

    double rvv10_energy = 0;
    double rvv10_vtxc = 0;
    ModuleBase::matrix rvv10_vofk;
    elecstate::PotXC base(&pw, &rvv10_energy, &rvv10_vtxc, &rvv10_vofk);
    elecstate::PotRvv10 rvv10(&pw, &rvv10_energy, &rvv10_vtxc);
    XC_Functional::set_xc_type("GGA_X_RPW86+GGA_C_PBE");
    ModuleBase::matrix rvv10_potential(1, pw.nrxx);
    base.cal_v_eff(&charge, &cell, rvv10_potential);
    rvv10.cal_v_eff(&charge, &cell, rvv10_potential);
    ASSERT_GT(std::abs(rvv10_energy - pbe_energy), 1.e-4);

    XC_Functional::set_xc_type("PBE");
    ModuleBase::matrix restored(1, pw.nrxx);
    pbe.cal_v_eff(&charge, &cell, restored);
    EXPECT_NEAR(pbe_energy, reference_pbe_energy, 1.e-12);
    EXPECT_NEAR(pbe_vtxc, reference_pbe_vtxc, 1.e-12);
    for (int i = 0; i < pw.nrxx; ++i)
        EXPECT_NEAR(restored(0, i), pbe_potential(0, i), 1.e-13);
}

} // namespace
