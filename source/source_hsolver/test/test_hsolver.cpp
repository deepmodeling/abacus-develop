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
        hsolver::set_diagethr_ks("lcao", "ksdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        0.0);
    EXPECT_EQ(
        hsolver::set_diagethr_ks("pw", "sdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        0.0);
}

TEST(HSolverSetDiagethrKsTest, Nscf)
{
    // the default 1e-2 is tightened to 0.1 * min(1e-2, scf_thr / nelec)
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "nscf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        1.25e-8);
    // a threshold the user already tightened is kept
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "nscf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6),
        1.0e-6);
}

TEST(HSolverSetDiagethrKsTest, FirstIteration)
{
    // atomic starting charge keeps the loose 1e-2
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        1.0e-2);
    // a charge read from file is trusted, so the first diagonalization is strict
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "file", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        1.0e-5);
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "wfc", "double", 0, 1, 1.0e-3, 1.0e-2, 1.0e-2, 8.0, 1.0e-6),
        1.0e-5);
    // relax/md never go below pw_diag_thr_init
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "relax", "file", "double", 0, 1, 1.0e-3, 1.0e-3, 1.0e-2, 8.0, 1.0e-6),
        1.0e-3);
}

TEST(HSolverSetDiagethrKsTest, LaterIterations)
{
    // iteration 2 restarts from 1e-2 and is bounded by 0.1 * drho / nelec
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 2, 1.0e-3, 1.0e-2, 1.0e-7, 8.0, 1.0e-6),
        1.25e-5);
    // later iterations only ever tighten the threshold
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 3, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6),
        1.0e-6);
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "double", 0, 3, 1.0e-3, 1.0e-2, 1.0e-4, 8.0, 1.0e-6),
        1.25e-5);
}

TEST(HSolverSetDiagethrKsTest, SinglePrecisionLowerBound)
{
    EXPECT_DOUBLE_EQ(
        hsolver::set_diagethr_ks("pw", "ksdft", "scf", "atomic", "single", 0, 3, 1.0e-3, 1.0e-2, 1.0e-6, 8.0, 1.0e-6),
        0.5e-4);
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
    const std::string log_name = "test_hsolver_reset_diag_ethr.log";
    std::ofstream ofs(log_name.c_str());
    const double new_ethr = hsolver::reset_diag_ethr(ofs, "pw", "ksdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0);
    const double new_ethr_single = hsolver::reset_diag_ethr(ofs, "pw", "ksdft", "single", 1.0e-2, 1.0e-3, 1.0e-2, 8.0);
    ofs.close();

    EXPECT_DOUBLE_EQ(new_ethr, 1.25e-5);
    EXPECT_DOUBLE_EQ(new_ethr_single, 0.5e-4);

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
    EXPECT_EQ(hsolver::reset_diag_ethr(ofs, "lcao", "ksdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0), 0.0);
    EXPECT_EQ(hsolver::reset_diag_ethr(ofs, "pw", "sdft", "double", 1.0e-2, 1.0e-3, 1.0e-2, 8.0), 0.0);
}
