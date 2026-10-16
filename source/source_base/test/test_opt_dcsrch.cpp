#include "../opt_dcsrch.h"
#include "gtest/gtest.h"

TEST(OptDCsrch, ConvergenceAndRestart)
{
    ModuleBase::Opt_DCsrch search;
    search.set_paras();
    std::string task;
    for (int run = 0; run < 2; ++run)
    {
        double f = 1.0;
        double g = -2.0;
        double step = 1.0;
        task = "START";
        search.dcSrch(f, g, step, task);
        ASSERT_EQ(task, "FG");
        EXPECT_EQ(task.size(), 2u);
        f = (step - 1.0) * (step - 1.0);
        g = 2.0 * (step - 1.0);
        search.dcSrch(f, g, step, task);
        EXPECT_EQ(task, "CONVERGENCE");
        EXPECT_EQ(task.size(), 11u);
        EXPECT_DOUBLE_EQ(step, 1.0);
    }
}

TEST(OptDCsrch, ErrorReplacesStatus)
{
    ModuleBase::Opt_DCsrch search;
    search.set_paras();
    double f = 1.0;
    double g = 1.0;
    double step = 1.0;
    std::string task = "START";
    search.dcSrch(f, g, step, task);
    EXPECT_EQ(task, "ERROR: INITIAL G .GE. ZERO");
    EXPECT_EQ(task.size(), std::string("ERROR: INITIAL G .GE. ZERO").size());

    task = "START";
    g = -2.0;
    search.dcSrch(f, g, step, task);
    EXPECT_EQ(task, "FG");
    EXPECT_EQ(task.size(), 2u);
}

TEST(OptDCsrch, WarningAtMaximumStep)
{
    ModuleBase::Opt_DCsrch search;
    search.set_paras(1e-4, 0.2, 1e-12, 0.0, 0.1);
    double f = 1.0;
    double g = -2.0;
    double step = 0.1;
    std::string task = "START";
    search.dcSrch(f, g, step, task);
    ASSERT_EQ(task, "FG");
    f = (step - 1.0) * (step - 1.0);
    g = 2.0 * (step - 1.0);
    search.dcSrch(f, g, step, task);
    EXPECT_EQ(task, "WARNING: STP = STPMAX");
    EXPECT_DOUBLE_EQ(step, 0.1);
}
