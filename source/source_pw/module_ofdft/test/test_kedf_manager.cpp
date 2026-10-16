#include "../kedf_manager.h"

#include "source_base/constants.h"
#include "source_io/module_parameter/input_parameter.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

class TestKEDFManager : public ::testing::Test
{
protected:
    void SetUp() override
    {
        const double length = ModuleBase::TWO_PI;
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
#ifdef __MPI
        basis.initmpi(1, 0, MPI_COMM_SELF);
#endif
        basis.initgrids(length, lattice, 8, 8, 8);
        basis.initparameters(false, 4.0, 1, false);
        basis.setuptransform();
        basis.collect_local_pw();
        volume = length * length * length;
        dV = volume / basis.nxyz;
        phi.resize(basis.nrxx);
        rho.resize(basis.nrxx);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const int iz = ir % basis.nplane + basis.startz_current;
            const double phase = ModuleBase::TWO_PI * iz / basis.nz;
            phi[ir] = 2.0 + amplitude * std::cos(phase);
            rho[ir] = phi[ir] * phi[ir];
            nelec += rho[ir] * dV;
        }
    }

    double analytic_tf_energy() const
    {
        double integral = 0.0;
        for (const double value : rho)
        {
            integral += std::pow(value, 5.0 / 3.0) * dV;
        }
        return tf_coefficient * integral;
    }

    double analytic_vw_energy() const
    {
        return volume * amplitude * amplitude / 2.0;
    }

    void check_wt(const Input_para& inp)
    {
        const double* density[] = {rho.data()};
        const double* wave[] = {phi.data()};
        double* mutable_wave[] = {phi.data()};
        KEDF_TF tf;
        KEDF_vW vw;
        KEDF_WT wt;
        tf.set_para(basis.nrxx, dV, inp.of_tf_weight);
        vw.set_para(dV, inp.of_vw_weight);
        wt.set_para(dV, inp.of_wt_alpha, inp.of_wt_beta, nelec,
                    inp.of_tf_weight, inp.of_vw_weight, inp.of_wt_rho0,
                    inp.of_hold_rho0, inp.of_read_kernel, inp.of_kernel_file, &basis);
        const double tf_energy = tf.get_energy(density);
        const double vw_energy = vw.get_energy(mutable_wave, &basis);
        const double wt_energy = wt.get_energy(density, &basis);
        ASSERT_GT(std::abs(wt_energy), 1e-8);

        // Compare public dispatch with independently initialized components.
        ModuleBase::matrix expected(1, basis.nrxx);
        tf.tf_potential(density, expected);
        wt.wt_potential(density, &basis, expected);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            expected(0, ir) *= 2.0 * phi[ir];
        }
        vw.vw_potential(wave, &basis, expected);

        manager.init(inp, &basis, dV, nelec);
        ModuleBase::matrix potential(1, basis.nrxx);
        manager.get_potential(density, wave, &basis, potential);
        const double total = tf_energy + vw_energy + wt_energy;
        EXPECT_NEAR(manager.get_energy(), total, 1e-9);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            EXPECT_NEAR(potential(0, ir), expected(0, ir), 1e-10);
        }

        std::vector<std::string> titles;
        std::vector<double> energies;
        manager.record_energy(titles, energies);
        const std::vector<std::string> expected_titles = {"TF KEDF", "vW KEDF", "WT KEDF"};
        EXPECT_EQ(titles, expected_titles);
        ASSERT_EQ(energies.size(), 3u);
        EXPECT_NEAR(energies[0], tf_energy, 1e-9);
        EXPECT_NEAR(energies[1], vw_energy, 1e-9);
        EXPECT_NEAR(energies[2], wt_energy, 1e-9);
    }

    ModulePW::PW_Basis basis{"cpu", "double"};
    KEDF_Manager manager;
    std::vector<double> phi;
    std::vector<double> rho;
    double volume = 0.0;
    double dV = 0.0;
    double nelec = 0.0;
    const double amplitude = 0.25;
    const double tf_coefficient = 5.742468000376382;
};

TEST_F(TestKEDFManager, InitTFOnly)
{
    Input_para inp;
    inp.of_kinetic = "tf";
    manager.init(inp, &basis, dV, nelec);
    const double* density[] = {rho.data()};
    const double* wave[] = {phi.data()};
    ModuleBase::matrix potential(1, basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        potential(0, ir) = 0.3;
    }
    manager.get_potential(density, wave, &basis, potential);
    const double energy = analytic_tf_energy();
    EXPECT_NEAR(manager.get_energy(), energy, 1e-9);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double tf_potential = 5.0 / 3.0 * tf_coefficient * std::pow(rho[ir], 2.0 / 3.0);
        const double expected = (0.3 + tf_potential) * 2.0 * phi[ir];
        EXPECT_NEAR(potential(0, ir), expected, 1e-10);
    }
    std::vector<std::string> titles;
    std::vector<double> energies;
    manager.record_energy(titles, energies);
    const std::vector<std::string> expected_titles = {"TF KEDF"};
    EXPECT_EQ(titles, expected_titles);
    ASSERT_EQ(energies.size(), 1u);
    EXPECT_NEAR(energies[0], energy, 1e-9);
}

TEST_F(TestKEDFManager, InitVWOnly)
{
    Input_para inp;
    inp.of_kinetic = "vw";
    manager.init(inp, &basis, dV, nelec);
    const double* density[] = {rho.data()};
    const double* wave[] = {phi.data()};
    ModuleBase::matrix potential(1, basis.nrxx);
    manager.get_potential(density, wave, &basis, potential);
    const double energy = analytic_vw_energy();
    EXPECT_NEAR(manager.get_energy(), energy, 1e-10);
    std::vector<std::string> titles;
    std::vector<double> energies;
    manager.record_energy(titles, energies);
    const std::vector<std::string> expected_titles = {"vW KEDF"};
    EXPECT_EQ(titles, expected_titles);
    ASSERT_EQ(energies.size(), 1u);
    EXPECT_NEAR(energies[0], energy, 1e-10);
}

TEST_F(TestKEDFManager, InitTFPlusCombinesTFAndVW)
{
    Input_para inp;
    inp.of_kinetic = "tf+";
    inp.of_vw_weight = 1.0 / 9.0;
    manager.init(inp, &basis, dV, nelec);
    const double* density[] = {rho.data()};
    const double* wave[] = {phi.data()};
    ModuleBase::matrix potential(1, basis.nrxx);
    manager.get_potential(density, wave, &basis, potential);
    const double tf_energy = analytic_tf_energy();
    const double vw_energy = analytic_vw_energy() * inp.of_vw_weight;
    const double expected = tf_energy + vw_energy;
    EXPECT_NEAR(manager.get_energy(), expected, 1e-9);
    std::vector<std::string> titles;
    std::vector<double> energies;
    manager.record_energy(titles, energies);
    const std::vector<std::string> expected_titles = {"TF KEDF", "vW KEDF"};
    EXPECT_EQ(titles, expected_titles);
    ASSERT_EQ(energies.size(), 2u);
    EXPECT_NEAR(energies[0], tf_energy, 1e-9);
    EXPECT_NEAR(energies[1], vw_energy, 1e-10);
}

TEST_F(TestKEDFManager, InitWTAndReinitializeWithChangedParameters)
{
    Input_para inp;
    inp.of_kinetic = "wt";
    inp.of_wt_alpha = 0.7;
    inp.of_wt_beta = 0.9;
    inp.of_vw_weight = 0.4;
    check_wt(inp);

    // Reuse the manager and verify set_para refreshes existing objects.
    inp.of_wt_alpha = 0.8;
    inp.of_wt_beta = 1.0;
    inp.of_vw_weight = 0.8;
    check_wt(inp);
}
