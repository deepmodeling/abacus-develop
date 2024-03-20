#include "gtest/gtest.h"
#include "../opt_CG.h"
#include "../opt_DCsrch.h"
#include "./opt_test_tools.h"

#define DOUBLETHRESHOLD 1e-8

class CG_test : public testing::Test
{
protected:
    ModuleBase::Opt_CG cg;
    ModuleBase::Opt_DCsrch ds;
    // LinearEqu le;
    // MinFunc mf;
    TestTools tools;
    int maxiter = 500;
    double step = 1.;
    double residual = 10.;
    double tol = 1e-5;
    int final_iter = 0;
    char *task = NULL;
    double *Ap = NULL;
    double *p = NULL;
    double *x = NULL;

    void SetUp()
    {
        cg.set_para(1.);
        cg.allocate(tools.nx);
        cg.init_b(tools.le.b);
        task = new char[60];
        Ap = new double[tools.nx];
        p = new double[tools.nx];
        x = new double[tools.nx];
    }

    void TearDown()
    {
        delete[] task;
        delete[] Ap;
        delete[] p;
        delete[] x;
    }

    void CG_Solve_LinearEq()
    {
        final_iter = 0;
        cg.refresh(0, tools.le.b);
        step = 1;
        residual = 10.;
        for (int i = 0; i < tools.nx; ++i)
        {
            x[i] = 0;
            p[i] = 0;
            Ap[i] = 0;
        }
        for (int iter = 0; iter < maxiter; ++iter)
        {
            if (residual < tol) 
            {
                final_iter = iter;
                break;
            }
            cg.next_direct(Ap, 0, p);
            tools.le.get_Ap(tools.le.A, p, Ap);
            int ifPD = 0;
            step = cg.step_length(Ap, p, ifPD);
            for (int i = 0; i < 3; ++i) x[i] += step * p[i]; 
            residual = cg.get_residual();
        }
    }

    void Solve(int cg_label, int func_label)
    {
        if (func_label==0)
        {
            cg.refresh(0, tools.le.b);
        }
        else
        {
            cg.refresh();
        }
        ds.set_paras(1e-4, 2e-1, 1e-12, 0.,12.);
        step = 1.;
        residual = 10.;
        final_iter = 0;
        for (int i = 0; i < tools.nx; ++i)
        {
            x[i] = 0;
            p[i] = 0;
            Ap[i] = 0;
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
            cg.next_direct(gradient, cg_label, p);
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

TEST_F(CG_test, Stand_Solve_LinearEq)
{
    CG_Solve_LinearEq();
    EXPECT_NEAR(x[0], 0.5, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], 1.6429086563584579739e-18, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 1.5, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 4);
    ASSERT_EQ(cg.get_iter(), 4);
}

TEST_F(CG_test, PR_Solve_LinearEq)
{
    Solve(1, 0);
    EXPECT_NEAR(x[0], 0.50000000000003430589, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], -3.4028335704761047964e-14, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 1.5000000000000166533, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 3);
    ASSERT_EQ(cg.get_iter(), 3);
}

TEST_F(CG_test, HZ_Solve_LinearEq)
{
    Solve(2, 0);
    EXPECT_NEAR(x[0], 0.49999999999999944489, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], -9.4368957093138305936e-16, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 1.5000000000000011102, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 3);
    ASSERT_EQ(cg.get_iter(), 3);
}

TEST_F(CG_test, PR_Min_Func)
{
    Solve(1, 1);
    EXPECT_NEAR(x[0], 4.0006805979150792396, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], 2.0713759992720870429, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 9.2871067233169171118, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 18);
    ASSERT_EQ(cg.get_iter(), 18);
}

TEST_F(CG_test, HZ_Min_Func)
{
    Solve(2, 1);
    EXPECT_NEAR(x[0], 4.0006825378033568086, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[1], 2.0691732100663737803, DOUBLETHRESHOLD);
    EXPECT_NEAR(x[2], 9.2780872787668311474, DOUBLETHRESHOLD);
    ASSERT_EQ(final_iter, 18);
    ASSERT_EQ(cg.get_iter(), 18);
}
// g++ -std=c++11 ../opt_CG.cpp ../opt_DCsrch.cpp ./CG_test.cpp ./test_tools.cpp  -lgtest -lpthread -lgtest_main -o test.exe