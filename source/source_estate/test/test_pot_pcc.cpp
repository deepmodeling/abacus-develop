#include "source_estate/module_pot/pot_pcc.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <sstream>

// These fixtures use public cell/charge storage without running INPUT setup.
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
Charge::Charge() {}
Charge::~Charge() {}

namespace
{
class PotPccTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        basis.nx = 2;
        basis.ny = 2;
        basis.nz = 2;
        basis.nplane = 2;
        basis.startz_current = 0;
        basis.nrxx = 8;
        basis.nxyz = 8;
        cell.lat0 = 10.0;
        cell.omega = 1000.0;
        cell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, 1.0);
        cell.ntype = 1;
        cell.nat = 1;
        atom.na = 1;
        atom.mass = 1.0;
        atom.ncpp.zv = 1.0;
        atom.tau = {ModuleBase::Vector3<double>(0.45, 0.45, 0.45)};
        cell.atoms = &atom;
        density.assign(8, 0.0);
        density[7] = 1.0 / 125.0;
        channels[0] = density.data();
        charge.rho = channels;
        charge.nspin = 1;
    }

    ModulePW::PW_Basis basis{"cpu", "double"};
    UnitCell cell;
    Atom atom;
    Charge charge;
    std::vector<double> density;
    double* channels[2] = {nullptr, nullptr};
};
} // namespace

TEST_F(PotPccTest, NeutralDipoleProducesConsistentEnergyPotentialAndForce)
{
    elecstate::PotPcc correction(&basis, 1.0);
    ModuleBase::matrix potential(1, 8);
    correction.cal_v_eff(&charge, &cell, potential);
    const double expected_energy = 4.0 * ModuleBase::PI * 0.75 / 3000.0;
    EXPECT_NEAR(correction.get_energy(), expected_energy, 1.0e-14);
    EXPECT_NEAR(potential(0, 7), expected_energy, 1.0e-14);
    const std::vector<double>& stored = correction.electron_potential();
    ASSERT_EQ(stored.size(), 8u);
    EXPECT_DOUBLE_EQ(stored[7], potential(0, 7));
    ModuleBase::matrix force(1, 3);
    correction.add_force(cell, force);
    const double expected_force = 8.0 * ModuleBase::PI * 0.5 / 3000.0;
    EXPECT_NEAR(force(0, 0), expected_force, 1.0e-14);
    EXPECT_NEAR(force(0, 1), expected_force, 1.0e-14);
    EXPECT_NEAR(force(0, 2), expected_force, 1.0e-14);
}

TEST_F(PotPccTest, RecomputesDensityKeepsInstancesIndependentAndSharesSpinChannels)
{
    elecstate::PotPcc first(&basis, 1.0);
    elecstate::PotPcc second(&basis, 1.0);
    ModuleBase::matrix potential(1, 8);
    first.cal_v_eff(&charge, &cell, potential);
    const double neutral_energy = first.get_energy();
    density[7] = 2.0 / 125.0;
    potential.zero_out();
    second.cal_v_eff(&charge, &cell, potential);
    const double charged_energy = second.get_energy();
    const double expected_charged = 2.837297479480619 / 10.0 + 4.0 * ModuleBase::PI * 1.5 / 3000.0;
    EXPECT_NEAR(charged_energy, expected_charged, 1.0e-14);
    EXPECT_DOUBLE_EQ(first.get_energy(), neutral_energy);
    potential.zero_out();
    first.cal_v_eff(&charge, &cell, potential);
    EXPECT_NEAR(first.get_energy(), charged_energy, 1.0e-14);

    // Both spin channels get the same correction on top of their potentials.
    std::vector<double> down(8, 0.0);
    density[7] = 0.3 / 125.0;
    down[7] = 0.7 / 125.0;
    channels[1] = down.data();
    charge.nspin = 2;
    elecstate::PotPcc correction(&basis, 1.0);
    ModuleBase::matrix spin_potential(2, 8);
    spin_potential(0, 7) = 2.0;
    spin_potential(1, 7) = 3.0;
    correction.cal_v_eff(&charge, &cell, spin_potential);
    const double energy = correction.get_energy();
    const double first_expected = 2.0 + energy;
    const double second_expected = 3.0 + energy;
    EXPECT_NEAR(spin_potential(0, 7), first_expected, 1.0e-14);
    EXPECT_NEAR(spin_potential(1, 7), second_expected, 1.0e-14);
}

TEST_F(PotPccTest, ForceMatchesFixedDensityEnergyDerivativeAfterMovingIons)
{
    elecstate::PotPcc correction(&basis, 1.0);
    ModuleBase::matrix potential(1, 8);
    correction.cal_v_eff(&charge, &cell, potential);
    ModuleBase::matrix force(1, 3);
    correction.add_force(cell, force);
    const double step_bohr = 1.0e-4;
    const double step_lat0 = step_bohr / cell.lat0;
    atom.tau[0].x += step_lat0;
    potential.zero_out();
    correction.cal_v_eff(&charge, &cell, potential);
    const double plus_energy = correction.get_energy();
    atom.tau[0].x -= 2.0 * step_lat0;
    potential.zero_out();
    correction.cal_v_eff(&charge, &cell, potential);
    const double minus_energy = correction.get_energy();
    const double numerical_force = -(plus_energy - minus_energy) / (2.0 * step_bohr);
    EXPECT_NEAR(force(0, 0), numerical_force, 1.0e-10);
}

TEST_F(PotPccTest, SlabAxisControlsPotentialEnergyForceAndIgnoresInPlaneCoordinates)
{
    const double expected_energy = 4.0 * ModuleBase::PI * 0.25 / 1000.0;
    const double expected_force = 8.0 * ModuleBase::PI * 0.5 / 1000.0;
    for (int axis = 0; axis < 3; ++axis)
    {
        elecstate::PotPcc correction(&basis, elecstate::PotPcc::Dimension::slab, axis, 1.0);
        ModuleBase::matrix potential(1, 8);
        correction.cal_v_eff(&charge, &cell, potential);
        EXPECT_NEAR(correction.get_energy(), expected_energy, 1.0e-14);
        EXPECT_NEAR(potential(0, 7), expected_energy, 1.0e-14);
        ModuleBase::matrix force(1, 3);
        correction.add_force(cell, force);
        for (int component = 0; component < 3; ++component)
        {
            const double expected = (component == axis) ? expected_force : 0.0;
            EXPECT_NEAR(force(0, component), expected, 1.0e-14);
        }
    }
}

TEST_F(PotPccTest, RejectsUnsupportedCellsKpointsAndSymmetry)
{
    const std::vector<ModuleBase::Vector3<double>> points = {ModuleBase::Vector3<double>(0.25, 0.0, 0.0)};
    elecstate::PotPcc::validate_kpoints(points, 1, 2);
    EXPECT_EXIT(elecstate::PotPcc::validate_kpoints(points, 1, 0), testing::ExitedWithCode(1), "");
    elecstate::PotPcc correction(&basis, elecstate::PotPcc::Dimension::slab, 2, 1.0);
    cell.symm.ptrans = {ModuleBase::Vector3<double>(0.5, 0.0, 0.0)};
    ModuleBase::matrix potential(1, 8);
    correction.cal_v_eff(&charge, &cell, potential);
    cell.symm.ptrans[0].z = 0.5;
    EXPECT_EXIT(correction.cal_v_eff(&charge, &cell, potential), testing::ExitedWithCode(1), "");
    cell.symm.ptrans.clear();
    cell.symm.nrotk = 1;
    cell.symm.gmatrix[0] = ModuleBase::Matrix3(0.0, 0.0, 1.0,
                                             0.0, 1.0, 0.0,
                                             1.0, 0.0, 0.0);
    EXPECT_EXIT(correction.cal_v_eff(&charge, &cell, potential), testing::ExitedWithCode(1), "");
    cell.symm.nrotk = 0;
    cell.latvec.e22 = 2.0;
    elecstate::PotPcc molecule(&basis, 1.0);
    EXPECT_EXIT(molecule.cal_v_eff(&charge, &cell, potential), testing::ExitedWithCode(1), "");
    cell.latvec.e22 = 1.0;
    cell.latvec.e31 = 0.1;
    elecstate::PotPcc tilted_slab(&basis, elecstate::PotPcc::Dimension::slab, 2, 1.0);
    EXPECT_EXIT(tilted_slab.cal_v_eff(&charge, &cell, potential), testing::ExitedWithCode(1), "");
}

TEST_F(PotPccTest, ReportsGridChargeMismatchAndUsesGridCharge)
{
    const std::string log_name = "pot_pcc_charge_warning.log";
    std::ofstream& warning_log = GlobalV::ofs_warning;
    warning_log.open(log_name.c_str());
    ModuleBase::matrix matched_potential(1, 8);
    elecstate::PotPcc matched(&basis, 1.0);
    matched.cal_v_eff(&charge, &cell, matched_potential);
    warning_log.flush();
    std::ifstream matched_log(log_name.c_str());
    std::stringstream matched_text;
    matched_text << matched_log.rdbuf();
    EXPECT_EQ(matched_text.str().find("PCC net charge"), std::string::npos);

    ModuleBase::matrix mismatched_potential(1, 8);
    elecstate::PotPcc mismatched(&basis, 2.0);
    mismatched.cal_v_eff(&charge, &cell, mismatched_potential);
    warning_log.close();
    std::ifstream mismatched_log(log_name.c_str());
    std::stringstream mismatched_text;
    mismatched_text << mismatched_log.rdbuf();
    EXPECT_NE(mismatched_text.str().find("PCC net charge differs from the expected value"), std::string::npos);
    EXPECT_DOUBLE_EQ(mismatched.get_energy(), matched.get_energy());
    std::remove(log_name.c_str());
}
