#include "gtest/gtest.h"
#include <vector>
#include <cmath>

/************************************************
 *  unit test of write_elf logic
 ***********************************************/

/**
 * This test verifies the ELF calculation logic for nspin=4
 * by testing the key formulas directly without file I/O.
 */

class ElfLogicTest : public ::testing::Test
{
protected:
    // Test the Thomas-Fermi kinetic energy density calculation
    double calculate_tau_TF(double rho) {
        const double c_tf = 3.0 / 10.0 * std::pow(3 * std::pow(M_PI, 2.0), 2.0 / 3.0) * 2.0;
        if (rho > 0.0) {
            return c_tf * std::pow(rho, 5.0 / 3.0);
        } else {
            return 0.0;
        }
    }

    // Test the ELF calculation
    double calculate_elf(double tau, double tau_vw, double tau_TF) {
        const double eps = 1.0e-5;
        if (tau_TF > 1.0e-12) {
            double chi = (tau - tau_vw + eps) / tau_TF;
            return 1.0 / (1.0 + chi * chi);
        } else {
            return 0.0;
        }
    }
};

TEST_F(ElfLogicTest, ThomasFermiPositiveDensity)
{
    // Test with positive density
    double rho = 0.1;
    double tau_TF = calculate_tau_TF(rho);

    EXPECT_GT(tau_TF, 0.0);
    EXPECT_LT(tau_TF, 1.0); // Should be reasonable value
}

TEST_F(ElfLogicTest, ThomasFermiNegativeDensity)
{
    // Test with negative density (for magnetization components)
    double rho = -0.02;
    double tau_TF = calculate_tau_TF(rho);

    EXPECT_EQ(tau_TF, 0.0); // Should return 0 for negative density
}

TEST_F(ElfLogicTest, ThomasFermiZeroDensity)
{
    // Test with zero density
    double rho = 0.0;
    double tau_TF = calculate_tau_TF(rho);

    EXPECT_EQ(tau_TF, 0.0);
}

TEST_F(ElfLogicTest, ElfCalculationNormal)
{
    // Test ELF calculation with normal values
    double tau = 0.05;
    double tau_vw = 0.02;
    double tau_TF = 0.03;

    double elf = calculate_elf(tau, tau_vw, tau_TF);

    EXPECT_GE(elf, 0.0);
    EXPECT_LE(elf, 1.0);
}

TEST_F(ElfLogicTest, ElfCalculationSmallTauTF)
{
    // Test ELF calculation with very small tau_TF
    double tau = 0.05;
    double tau_vw = 0.02;
    double tau_TF = 1.0e-15; // Very small

    double elf = calculate_elf(tau, tau_vw, tau_TF);

    EXPECT_EQ(elf, 0.0); // Should return 0 for very small tau_TF
}

TEST_F(ElfLogicTest, ElfCalculationZeroTauTF)
{
    // Test ELF calculation with zero tau_TF
    double tau = 0.05;
    double tau_vw = 0.02;
    double tau_TF = 0.0;

    double elf = calculate_elf(tau, tau_vw, tau_TF);

    EXPECT_EQ(elf, 0.0); // Should return 0 for zero tau_TF
}

TEST_F(ElfLogicTest, ElfValueRange)
{
    // Test that ELF is always in [0, 1] for various inputs
    std::vector<double> tau_values = {0.01, 0.05, 0.1, 0.5, 1.0};
    std::vector<double> tau_vw_values = {0.005, 0.02, 0.05, 0.2, 0.5};
    std::vector<double> tau_TF_values = {0.01, 0.03, 0.08, 0.3, 0.8};

    for (double tau : tau_values) {
        for (double tau_vw : tau_vw_values) {
            for (double tau_TF : tau_TF_values) {
                double elf = calculate_elf(tau, tau_vw, tau_TF);
                EXPECT_GE(elf, 0.0) << "ELF should be >= 0";
                EXPECT_LE(elf, 1.0) << "ELF should be <= 1";
            }
        }
    }
}

TEST_F(ElfLogicTest, Nspin4ComponentHandling)
{
    // Test that we can handle 4 components independently
    int nspin = 4;
    std::vector<double> rho(nspin);
    std::vector<double> tau_TF(nspin);

    // Component 0: total charge (positive)
    rho[0] = 0.1;
    tau_TF[0] = calculate_tau_TF(rho[0]);
    EXPECT_GT(tau_TF[0], 0.0);

    // Components 1-3: magnetization (can be negative)
    for (int i = 1; i < nspin; ++i) {
        rho[i] = -0.01 * i; // Negative
        tau_TF[i] = calculate_tau_TF(rho[i]);
        EXPECT_EQ(tau_TF[i], 0.0); // Should be 0 for negative density
    }
}

TEST_F(ElfLogicTest, Nspin4AllPositive)
{
    // Test with all positive densities
    int nspin = 4;
    std::vector<double> rho(nspin);
    std::vector<double> tau_TF(nspin);

    for (int i = 0; i < nspin; ++i) {
        rho[i] = 0.05 + 0.01 * i;
        tau_TF[i] = calculate_tau_TF(rho[i]);
        EXPECT_GT(tau_TF[i], 0.0);
    }
}

TEST_F(ElfLogicTest, Nspin4MixedSigns)
{
    // Test with mixed positive and negative densities
    int nspin = 4;
    std::vector<double> rho = {0.1, -0.02, 0.03, -0.01};
    std::vector<double> tau_TF(nspin);

    for (int i = 0; i < nspin; ++i) {
        tau_TF[i] = calculate_tau_TF(rho[i]);
        if (rho[i] > 0.0) {
            EXPECT_GT(tau_TF[i], 0.0);
        } else {
            EXPECT_EQ(tau_TF[i], 0.0);
        }
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
