#include "../sto_dos_quality.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

/************************************************
 *  unit test of sto_dos_quality.cpp
 *
 *  Regression tests for the two quality guards of the SDFT stochastic DOS:
 *
 *  (a) the Chebyshev energy window must cover the deterministic eigenvalues,
 *      the plane-wave spectrum and the requested output grid - otherwise a
 *      positive lower bound hides every occupied state, and |t| > 1 makes
 *      T_n(t) grow exponentially (the historical 1e26--1e64 DOS values);
 *
 *  (b) DOS parameters that are silently useless at high temperature (grid
 *      spacing vs Gaussian broadening, Chebyshev order vs range/sigma) and a
 *      DOS curve that does not integrate to the number of plane-wave states,
 *      goes negative, or was evaluated outside the window.
 *
 *  All assertions encode the contract only, never an implementation detail, so
 *  any correct implementation passes.  The helpers only report warnings, hence
 *  "is a warning produced" / "is a healthy setup silent" is what is checked.
 ***********************************************/

namespace
{
const double ecut = 900.0;       // ecutwfc of the production run [Ry]
const double bad_emin = 112.433; // positive automatic lower bound of that run [Ry]
const double bad_emax = 810.199; // upper bound of the same run [Ry]
const double ks_min = -32.8;     // lowest true KS eigenvalue [Ry]
const double ks_max = 9.65;      // highest true KS eigenvalue [Ry]

std::vector<double> ks_eigs()
{
    return {ks_min, -12.0, -3.2, 0.0, 4.4, ks_max};
}

bool has_warning_containing(const std::vector<std::string>& warnings, const std::string& needle)
{
    for (const std::string& w: warnings)
    {
        if (w.find(needle) != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

/// Contract for a repaired window: finite, ordered, covering the eigenvalues
/// and the plane-wave spectrum top.
void expect_healthy(const ModuleIO::CheRange& r, const std::vector<double>& det_eigs, double ecut_in)
{
    EXPECT_TRUE(std::isfinite(r.emin));
    EXPECT_TRUE(std::isfinite(r.emax));
    EXPECT_GT(r.emax - r.emin, 1e-6);
    if (ecut_in > 0.0)
    {
        EXPECT_GE(r.emax, 1.05 * ecut_in);
    }
    for (double e: det_eigs)
    {
        EXPECT_LT(r.emin, e);
        EXPECT_GT(r.emax, e);
    }
}
} // namespace

class TestStoDosQuality : public ::testing::Test
{
};

// ---------------------------------------------------------------------------
// (a) Chebyshev-window guard
// ---------------------------------------------------------------------------

// 1) healthy input must be returned bit for bit unchanged (no regression for
//    already working runs)
TEST_F(TestStoDosQuality, healthy_range_is_returned_unchanged)
{
    const double emin = -40.0; // below min(ks_eigs) = -32.8
    const double emax = 950.0; // above 1.05*ecutwfc = 945
    const ModuleIO::CheRange r = ModuleIO::sanitize_che_range(emin, emax, ecut, ks_eigs());
    EXPECT_EQ(r.emin, emin);
    EXPECT_EQ(r.emax, emax);
    EXPECT_FALSE(r.emin_repaired);
    EXPECT_FALSE(r.emax_repaired);
    EXPECT_FALSE(r.degenerate);
}

// 2) key regression: a positive lower bound must be repaired, otherwise every
//    bound state falls outside the window
TEST_F(TestStoDosQuality, positive_emin_is_repaired)
{
    const ModuleIO::CheRange r = ModuleIO::sanitize_che_range(bad_emin, bad_emax, ecut, {});
    EXPECT_TRUE(r.emin_repaired);
    EXPECT_LT(r.emin, 0.0);         // the original keeps +112.4 -> FAIL
    EXPECT_GE(r.emax, 1.05 * ecut); // 810.2 -> 945
    expect_healthy(r, {}, ecut);
}

// 3) Emax below the plane-wave spectrum top must be raised to 1.05*ecutwfc
TEST_F(TestStoDosQuality, emax_floor_from_ecut)
{
    const ModuleIO::CheRange r = ModuleIO::sanitize_che_range(-10.0, 700.0, ecut, {});
    EXPECT_GE(r.emax, 1.05 * ecut);
    EXPECT_TRUE(r.emax_repaired);
    // no floor is imposed when ecutwfc is unknown (0)
    const ModuleIO::CheRange r0 = ModuleIO::sanitize_che_range(-10.0, 700.0, 0.0, {});
    EXPECT_EQ(r0.emax, 700.0);
}

// 4) deterministic (KS) eigenvalues must be covered
TEST_F(TestStoDosQuality, det_eigs_are_covered)
{
    const ModuleIO::CheRange r1 = ModuleIO::sanitize_che_range(5.0, 950.0, ecut, ks_eigs());
    EXPECT_LT(r1.emin, ks_min);
    expect_healthy(r1, ks_eigs(), ecut);

    const ModuleIO::CheRange r2 = ModuleIO::sanitize_che_range(bad_emin, bad_emax, ecut, ks_eigs());
    EXPECT_LT(r2.emin, ks_min);
    expect_healthy(r2, ks_eigs(), ecut);

    const ModuleIO::CheRange r3 = ModuleIO::sanitize_che_range(-10.0, 8.0, ecut, {ks_max});
    EXPECT_GT(r3.emax, ks_max);
}

// 5) invalid input (NaN / inverted / zero width) must become a usable window
TEST_F(TestStoDosQuality, invalid_input_is_sanitized)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<std::pair<double, double>> bad = {{950.0, -10.0}, {nan, 950.0}, {-10.0, nan}, {0.0, 0.0}};
    for (const std::pair<double, double>& p: bad)
    {
        const ModuleIO::CheRange r = ModuleIO::sanitize_che_range(p.first, p.second, ecut, {});
        EXPECT_TRUE(r.degenerate);
        expect_healthy(r, {}, ecut);
    }
}

// 6) idempotent: a repaired window must be stable, otherwise the SCF drifts
TEST_F(TestStoDosQuality, idempotent)
{
    const ModuleIO::CheRange r1 = ModuleIO::sanitize_che_range(bad_emin, bad_emax, ecut, ks_eigs());
    const ModuleIO::CheRange r2 = ModuleIO::sanitize_che_range(r1.emin, r1.emax, ecut, ks_eigs());
    EXPECT_EQ(r1.emin, r2.emin);
    EXPECT_EQ(r1.emax, r2.emax);
    EXPECT_FALSE(r2.emin_repaired);
    EXPECT_FALSE(r2.emax_repaired);
}

// 7) the DOS output window must stay inside the Chebyshev window so that
//    |t| <= 1 holds by construction (dos_scale used to break this)
TEST_F(TestStoDosQuality, healthy_parameters_are_silent)
{
    // range 2260 eV (typical high-T case), sigma 20 eV, de 1 eV, nche 600
    EXPECT_TRUE(ModuleIO::check_dos_parameters(2260.0, 20.0, 1.0, 600).empty());
    EXPECT_TRUE(ModuleIO::check_dos_parameters(40.0, 0.1, 0.05, 2000).empty());
}

// 9) grid vs broadening: de << sigma and de >> sigma both warn
TEST_F(TestStoDosQuality, grid_density_versus_sigma)
{
    const std::vector<std::string> fine = ModuleIO::check_dos_parameters(200.0, 1.0, 0.001, 500);
    EXPECT_TRUE(has_warning_containing(fine, "much smaller than dos_sigma"));

    const std::vector<std::string> coarse = ModuleIO::check_dos_parameters(200.0, 0.1, 5.0, 500);
    EXPECT_TRUE(has_warning_containing(coarse, "much larger than dos_sigma"));
}

// 10) dos_sigma <= 0 and dos_edelta_ev <= 0 are reported immediately
TEST_F(TestStoDosQuality, non_positive_sigma_is_reported)
{
    const std::vector<std::string> zero = ModuleIO::check_dos_parameters(100.0, 0.0, 1.0, 500);
    ASSERT_FALSE(zero.empty());
    EXPECT_TRUE(has_warning_containing(zero, "not positive"));
    EXPECT_TRUE(has_warning_containing(ModuleIO::check_dos_parameters(100.0, -1.0, 1.0, 500), "not positive"));

    const std::vector<std::string> bad_de = ModuleIO::check_dos_parameters(100.0, 1.0, 0.0, 500);
    ASSERT_FALSE(bad_de.empty());
    EXPECT_TRUE(has_warning_containing(bad_de, "dos_edelta_ev must be positive"));
}

// 11) insufficient dos_nche must report the estimated order needed
TEST_F(TestStoDosQuality, insufficient_nche_is_reported)
{
    const std::vector<std::string> warn = ModuleIO::check_dos_parameters(2260.0, 0.07, 0.05, 600);
    ASSERT_FALSE(warn.empty());
    EXPECT_TRUE(has_warning_containing(warn, "narrower than the Chebyshev expansion can resolve"));
    EXPECT_TRUE(has_warning_containing(warn, "estimated order needed"));
}

// 12) an extremely wide range relative to sigma is reported

TEST_F(TestStoDosQuality, result_checks_are_silent_for_a_healthy_curve)
{
    // |t| <= 1 and a positive curve: nothing to report.
    EXPECT_TRUE(ModuleIO::check_dos_result(0.5, 1.0).empty());
}

TEST_F(TestStoDosQuality, result_check_reports_energies_outside_the_interval)
{
    // One output energy at |t| = 1.5: T_n(t) diverges there.
    const std::vector<std::string> w = ModuleIO::check_dos_result(0.5, 1.5);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("outside the Chebyshev interval"), std::string::npos);
}

TEST_F(TestStoDosQuality, result_check_reports_a_negative_curve)
{
    const std::vector<std::string> w = ModuleIO::check_dos_result(-3.0, 0.9);
    ASSERT_EQ(w.size(), 1u);
    EXPECT_NE(w[0].find("negative"), std::string::npos);
}
