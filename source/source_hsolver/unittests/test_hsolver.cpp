#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include "source_hsolver/hsolver.h"

/************************************************
 *  unit test of the free functions in hsolver.h
 ***********************************************/

/**
 * Tested functions:
 *  - hsolver::cal_hsolve_error
 *      - pw + ksdft: diag_ethr * max(1, nelec)
 *      - any other basis/esolver: 0
 *  - hsolver::set_diagethr_ks
 *      - non pw-ksdft returns 0
 *      - nscf, first scf iteration (by init_chg and calculation), later iterations
 *      - single precision lower bound
 *      - optional SCF/relax upper bound, including every electronic-cycle reset
 *  - hsolver::set_diagethr_sdft
 *      - non pw-sdft returns 0
 *      - nscf, first iteration of the first/later ion step, later iterations
 *  - hsolver::reset_diag_ethr
 *      - pw + ksdft: 0.1 * drho / nelec, with the single precision lower bound
 *      - any other basis/esolver: 0
 */

TEST(HSolverCalHsolveErrorTest, PwKsdftScalesWithNelec)
{
    EXPECT_DOUBLE_EQ(hsolver::cal_hsolve_error("pw", "ksdft", 1.0e-3, 8.0), 8.0e-3);
    // nelec below 1 is clamped to 1
    EXPECT_DOUBLE_EQ(hsolver::cal_hsolve_error("pw", "ksdft", 1.0e-3, 0.5), 1.0e-3);
}

TEST(HSolverCalHsolveErrorTest, OtherCasesAreZero)
{
    EXPECT_EQ(hsolver::cal_hsolve_error("lcao", "ksdft", 1.0e-3, 8.0), 0.0);
    EXPECT_EQ(hsolver::cal_hsolve_error("lcao_in_pw", "ksdft", 1.0e-3, 8.0), 0.0);
    EXPECT_EQ(hsolver::cal_hsolve_error("pw", "sdft", 1.0e-3, 8.0), 0.0);
}

TEST(HSolverSetDiagethrKsTest, NonPwKsdftIsZero)
{
    EXPECT_EQ(
        hsolver::set_diagethr_ks("lcao", "ksdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        0.0);
    EXPECT_EQ(
        hsolver::set_diagethr_ks("pw", "sdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        0.0);
}

TEST(HSolverSetDiagethrKsTest, Nscf)
{
    // the default 1e-2 is tightened to 0.1 * min(1e-2, scf_thr / nelec)
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "nscf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        1.25e-8);
    // a threshold the user already tightened is kept
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "nscf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6, 0.0),
        1.0e-6);
}

TEST(HSolverSetDiagethrKsTest, FirstIteration)
{
    // atomic starting charge keeps the loose 1e-2
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        1.0e-2);
    // a charge read from file is trusted, so the first diagonalization is strict
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "file", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        1.0e-5);
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "wfc", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 0.0),
        1.0e-5);
    // relax/md never go below pw_diag_thr_init
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "relax", "file", "double", 0, 1, 1.0e-3, 1.0e-3, 1.0e-2, 8.0, 1.0e-6, 0.0),
        1.0e-3);
}

TEST(HSolverSetDiagethrKsTest, LaterIterations)
{
    // iteration 2 restarts from 1e-2 and is bounded by 0.1 * drho / nelec
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 2, 1.0e-3, 1.0e-2, 1.0e-7, 8.0, 1.0e-6, 0.0),
        1.25e-5);
    // later iterations only ever tighten the threshold
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 3, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6, 0.0),
        1.0e-6);
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 3, 1.0e-3, 1.0e-2, 1.0e-4, 8.0, 1.0e-6, 0.0),
        1.25e-5);
}

TEST(HSolverSetDiagethrKsTest, SinglePrecisionLowerBound)
{
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "single", 0, 3, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6, 0.0),
        0.5e-4);
}

TEST(HSolverSetDiagethrKsTest, ScfUpperBoundAppliesFromFirstIteration)
{
    for (const int iter : {1, 2, 3})
    {
        const double threshold = hsolver::set_diagethr_ks(
            "pw", "ksdft", "scf", "atomic", "double", 0, iter, 1.0, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
        EXPECT_DOUBLE_EQ(threshold, 1.0e-12);
    }
}

TEST(HSolverSetDiagethrKsTest, ScfUpperBoundKeepsTighterThresholds)
{
    const double initial = hsolver::set_diagethr_ks(
        "pw", "ksdft", "scf", "atomic", "double", 0, 1, 1.0, 1.0e-14, 1.0e-14, 8.0, 1.0e-6, 1.0e-12);
    const double adaptive = hsolver::set_diagethr_ks(
        "pw", "ksdft", "scf", "atomic", "double", 0, 2, 8.0e-14, 1.0e-2, 1.0e-12, 8.0, 1.0e-6, 1.0e-12);
    const double previous = hsolver::set_diagethr_ks(
        "pw", "ksdft", "scf", "atomic", "double", 0, 3, 1.0, 1.0e-2, 1.0e-14, 8.0, 1.0e-6, 1.0e-12);
    EXPECT_DOUBLE_EQ(initial, 1.0e-14);
    EXPECT_DOUBLE_EQ(adaptive, 1.0e-15);
    EXPECT_DOUBLE_EQ(previous, 1.0e-14);
}

TEST(HSolverSetDiagethrKsTest, ScfUpperBoundPreservesSinglePrecisionFloor)
{
    const double capped = hsolver::set_diagethr_ks(
        "pw", "ksdft", "scf", "atomic", "single", 0, 2, 1.0, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-4);
    const double floor = hsolver::set_diagethr_ks(
        "pw", "ksdft", "scf", "atomic", "single", 0, 3, 1.0e-12, 1.0e-2, 1.0e-8, 8.0, 1.0e-6, 5.0e-5);
    EXPECT_DOUBLE_EQ(capped, 1.0e-4);
    EXPECT_DOUBLE_EQ(floor, 5.0e-5);
}

TEST(HSolverSetDiagethrKsTest, RelaxUpperBoundSurvivesEveryIonicCycle)
{
    // Carry the preceding electronic threshold into each new ionic cycle.
    // Both the first-iteration initial floor and iteration-2 reset can loosen it.
    for (const std::string& init_chg : {"atomic", "file"})
    {
        double previous = 1.0e-2;
        for (const int ionic_step : {0, 1, 2})
        {
            for (const int iteration : {1, 2, 3})
            {
                SCOPED_TRACE(::testing::Message() << init_chg << " ionic=" << ionic_step
                                                  << " electronic=" << iteration);
                const double threshold = hsolver::set_diagethr_ks(
                    "pw", "ksdft", "relax", init_chg, "double", ionic_step, iteration,
                    1.0, 1.0e-2, previous, 8.0, 1.0e-6, 1.0e-12);
                EXPECT_DOUBLE_EQ(threshold, 1.0e-12);
                previous = threshold;
            }
        }
    }
}

TEST(HSolverSetDiagethrKsTest, RelaxZeroUpperBoundPreservesLegacyIonicResets)
{
    const double residuals[] = {1.0, 8.0e-3, 8.0e-4};
    const double expected[] = {1.0e-3, 1.0e-4, 1.0e-5};
    double previous = 1.0e-2;
    for (const int ionic_step : {0, 1, 2})
    {
        for (int index = 0; index < 3; ++index)
        {
            const int iteration = index + 1;
            SCOPED_TRACE(::testing::Message() << "ionic=" << ionic_step << " electronic=" << iteration);
            const double threshold = hsolver::set_diagethr_ks(
                "pw", "ksdft", "relax", "file", "double", ionic_step, iteration,
                residuals[index], 1.0e-3, previous, 8.0, 1.0e-6, 0.0);
            EXPECT_DOUBLE_EQ(threshold, expected[index]);
            previous = threshold;
        }
    }
}

TEST(HSolverSetDiagethrKsTest, RelaxUpperBoundKeepsTighterThresholds)
{
    for (const int ionic_step : {0, 1, 2})
    {
        SCOPED_TRACE(ionic_step);
        const double initial = hsolver::set_diagethr_ks(
            "pw", "ksdft", "relax", "atomic", "double", ionic_step, 1,
            1.0, 1.0e-14, 1.0e-14, 8.0, 1.0e-6, 1.0e-12);
        const double adaptive = hsolver::set_diagethr_ks(
            "pw", "ksdft", "relax", "atomic", "double", ionic_step, 2,
            8.0e-14, 1.0e-2, 1.0e-12, 8.0, 1.0e-6, 1.0e-12);
        const double previous = hsolver::set_diagethr_ks(
            "pw", "ksdft", "relax", "atomic", "double", ionic_step, 3,
            1.0, 1.0e-2, 1.0e-14, 8.0, 1.0e-6, 1.0e-12);
        EXPECT_DOUBLE_EQ(initial, 1.0e-14);
        EXPECT_DOUBLE_EQ(adaptive, 1.0e-15);
        EXPECT_DOUBLE_EQ(previous, 1.0e-14);
    }
}

TEST(HSolverSetDiagethrKsTest, RelaxUpperBoundPreservesSinglePrecisionFloor)
{
    const double capped = hsolver::set_diagethr_ks(
        "pw", "ksdft", "relax", "atomic", "single", 1, 2, 1.0, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-4);
    const double floor = hsolver::set_diagethr_ks(
        "pw", "ksdft", "relax", "atomic", "single", 1, 3, 1.0e-12, 1.0e-2, 1.0e-8, 8.0, 1.0e-6, 5.0e-5);
    EXPECT_DOUBLE_EQ(capped, 1.0e-4);
    EXPECT_DOUBLE_EQ(floor, 5.0e-5);
}

TEST(HSolverSetDiagethrKsTest, UpperBoundDoesNotChangeOtherPolicies)
{
    const double nscf = hsolver::set_diagethr_ks(
        "pw", "ksdft", "nscf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
    const double cell_relax = hsolver::set_diagethr_ks(
        "pw", "ksdft", "cell-relax", "file", "double", 0, 1, 1.0e-3, 1.0e-3, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
    const double md = hsolver::set_diagethr_ks(
        "pw", "ksdft", "md", "file", "double", 0, 1, 1.0e-3, 1.0e-3, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
    const double lcao = hsolver::set_diagethr_ks(
        "lcao", "ksdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
    const double sdft = hsolver::set_diagethr_ks(
        "pw", "sdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6, 1.0e-12);
    EXPECT_DOUBLE_EQ(nscf, 1.25e-8);
    EXPECT_DOUBLE_EQ(cell_relax, 1.0e-3);
    EXPECT_DOUBLE_EQ(md, 1.0e-3);
    EXPECT_DOUBLE_EQ(lcao, 0.0);
    EXPECT_DOUBLE_EQ(sdft, 0.0);
}

TEST(HSolverSetDiagethrSdftTest, NonPwSdftIsZero)
{
    EXPECT_EQ(
        hsolver::set_diagethr_sdft("pw", "ksdft", "scf", "atomic", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 4, 2.0, 8.0, 1.0e-6),
        0.0);
    EXPECT_EQ(
        hsolver::set_diagethr_sdft("lcao", "sdft", "scf", "atomic", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 4, 2.0, 8.0, 1.0e-6),
        0.0);
}

TEST(HSolverSetDiagethrSdftTest, Nscf)
{
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "nscf", "atomic", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 4, 2.0, 8.0, 1.0e-6),
        1.25e-8);
}

TEST(HSolverSetDiagethrSdftTest, FirstIteration)
{
    // first ion step: file charge gives 1e-5, then bounded below by pw_diag_thr_init
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "scf", "file", 0, 1, 1.0e-3, 1.0e-2, 1.0e-4, 4, 2.0, 8.0, 1.0e-6),
        1.0e-2);
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "scf", "atomic", 0, 1, 1.0e-3, 1.0e-6, 1.0e-4, 4, 2.0, 8.0, 1.0e-6),
        1.0e-4);
    // later ion steps are bounded below by 1e-5
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "scf", "atomic", 1, 1, 1.0e-3, 1.0e-2, 1.0e-7, 4, 2.0, 8.0, 1.0e-6),
        1.0e-5);
}

TEST(HSolverSetDiagethrSdftTest, LaterIterations)
{
    // bounded by 0.1 * drho / KS electrons
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "scf", "atomic", 0, 2, 1.0e-3, 1.0e-2, 1.0e-2, 4, 2.0, 8.0, 1.0e-6),
        5.0e-5);
    // pure stochastic (no KS bands) does not diagonalize
    EXPECT_EQ(
        hsolver::set_diagethr_sdft("pw", "sdft", "scf", "atomic", 0, 2, 1.0e-3, 1.0e-2, 1.0e-2, 0, 2.0, 8.0, 1.0e-6),
        0.0);
}

TEST(HSolverResetDiagEthrTest, PwKsdft)
{
    const std::string log_name = testing::TempDir() + "test_hsolver_reset_diag_ethr.log";
    std::ofstream ofs(log_name.c_str());
    const double new_ethr = hsolver::reset_diag_ethr(ofs, "pw", "ksdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0, 0.0);
    const double new_ethr_single = hsolver::reset_diag_ethr(ofs, "pw", "ksdft", "single", 1.0e-2, 1.0e-3, 1.0e-2, 8.0, 0.0);
    // A fractional electron count can otherwise make the reset exceed the bound.
    const double capped = hsolver::reset_diag_ethr(
        ofs, "pw", "ksdft", "double", 1.0e-8, 5.0e-9, 1.0e-8, 1.0e-2, 1.0e-8);
    const double tighter = hsolver::reset_diag_ethr(
        ofs, "pw", "ksdft", "double", 1.0e-2, 1.0e-12, 1.0e-8, 8.0, 1.0e-8);
    const double single = hsolver::reset_diag_ethr(
        ofs, "pw", "ksdft", "single", 1.0e-2, 1.0e-12, 1.0e-8, 8.0, 5.0e-5);
    ofs.close();

    EXPECT_DOUBLE_EQ(new_ethr, 1.25e-5);
    EXPECT_DOUBLE_EQ(new_ethr_single, 0.5e-4);
    EXPECT_DOUBLE_EQ(capped, 1.0e-8);
    EXPECT_DOUBLE_EQ(tighter, 1.25e-14);
    EXPECT_DOUBLE_EQ(single, 5.0e-5);

    std::ifstream ifs(log_name.c_str());
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    ifs.close();
    EXPECT_NE(buffer.str().find("Threshold on eigenvalues was too large"), std::string::npos);
    EXPECT_NE(buffer.str().find("New diag ethr"), std::string::npos);
    std::remove(log_name.c_str());
}

TEST(HSolverResetDiagEthrTest, OtherCasesAreZero)
{
    std::ofstream ofs;
    EXPECT_EQ(hsolver::reset_diag_ethr(ofs, "lcao", "ksdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0, 0.0), 0.0);
    EXPECT_EQ(hsolver::reset_diag_ethr(ofs, "pw", "sdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0, 0.0), 0.0);
}
