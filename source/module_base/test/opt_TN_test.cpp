#include "gtest/gtest.h"
#include "../opt_TN.hpp"
#include "../opt_DCsrch.h"
#include "./opt_test_tools.h"

#define DOUBLETHRESHOLD 1e-5

class TN_test : public testing::Test
{
protected:
    ModuleBase::Opt_TN tn;
    ModuleBase::Opt_DCsrch ds;
    TestTools tools;
    int maxiter = 500;
    double step = 1.;
    double residual = 10.;
    double tol = 1e-5;
    int final_iter = 0;
    int flag = 0;
    char *task = NULL;
    double *p = NULL;
    double *x = NULL;

    void SetUp()
    {
        tn.set_para(1.);
        tn.allocate(tools.nx);
        task = new char[60];
        p = new double[tools.nx];
        x = new double[tools.nx];
    }

    void TearDown()
    {
        delete[] task;
        delete[] p;
        delete[] x;
    }

    void Solve(int func_label)
    {
        tn.refresh();
        ds.set_paras(1e-4, 2e-1, 1e-12, 0.,12.);
        step = 1.;
        residual = 10.;
        final_iter = 0;
        for (int i = 0; i < tools.nx; ++i)
        {
            x[i] = 0;
            p[i] = 0;
        }

        double f = 0;
        double g = 0;
        double *gradient = new double[3];
        double *temp_x = new double[3];
        ModuleBase::GlobalFunc::ZEROS(gradient, 3);
        ModuleBase::GlobalFunc::ZEROS(temp_x, 3);

        for (int iter = 0; iter < maxiter; ++iter)
        {
            tools.dfuncdx(x, gradient, func_label);
            residual = 0;
            for (int i = 0; i<3 ;++i) residual += gradient[i] * gradient[i];
            if (residual < tol) 
            {
                final_iter = iter;
                break;
            }
            if (func_label == 0)
            {
                tn.next_direct(x, gradient, flag, p, &(tools.le), &LinearEqu::dfuncdx);
            }
            else if (func_label == 1)
            {
                tn.next_direct(x, gradient, flag, p, &(tools.mf), &ModuleESolver::ESolver_OF::dfuncdx);
            }
            for (int i = 0; i < 3; ++i) temp_x[i] = x[i];
            task[0] = 'S'; task[1] = 'T'; task[2] = 'A'; task[3] = 'R'; task[4] = 'T';
            while (true)
            {
                f = tools.func(temp_x, func_label);
                g = tools.dfuncdstp(temp_x, p, func_label);
                ds.dcSrch(f, g, step, task);
                if (task[0] == 'F' && task[1] == 'G')
                {
                    for (int j = 0; j < 3; ++j) temp_x[j] = x[j] + step * p[j];
                    continue;
                }
                else if (task[0] == 'C' && task[1] == 'O')
                {
                    break;
                }
                else if (task[0] == 'W' && task[1] == 'A')
                {
                    break;
                } 
                else if (task[0] == 'E' && task[1] == 'R')
                {
                    break;
                }
            }
            for (int i = 0; i < 3; ++i) x[i] += step * p[i];
        }
        delete[] temp_x;
        delete[] gradient;
    }
};


TEST_F(TN_test, TN_Solve_LinearEq)
{
    Solve(0);
    EXPECT_NEAR(x[0], 0.50000000000003430589, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], -3.4028335704761047964e-14, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 1.5000000000000166533, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 1);
    ASSERT_EQ(tn.get_iter(), 1);
}

TEST_F(TN_test, TN_Min_Func)
{
    Solve(1);
    EXPECT_NEAR(x[0], 4.0049968540891525137, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], 2.1208751163987624722, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 9.4951527720891863993, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 6);
    ASSERT_EQ(tn.get_iter(), 6);
}