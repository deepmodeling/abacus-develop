#include "gtest/gtest.h"
#include "../broyden_mixing.h"
#include "../pulay_mixing.h"
#include "../plain_mixing.h"

#define DOUBLETHRESHOLD 1e-8

class Mixing_Test : public testing::Test
{
protected:
    const double mixing_beta = 0.8;
    const int mixing_ndim = 7;
    Base_Mixing::Mixing_Data xdata;
    Base_Mixing::Mixing *mixing;

    double thr = 1e-12;
    int niter = 0;
    int maxiter = 50;
    double* x_out = nullptr;

    void SetUp()
    {
        this->mixing = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    void TearDown()
    {
        delete this->mixing;
        delete[] this->x_out;
    }

    void solve(std::string method)
    {
        if (method == "broyden")
        {
            this->mixing = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
            this->xdata.resize(this->mixing_ndim, 3, sizeof(double));
        }
        else if (method == "pulay")
        {
            this->mixing = new Base_Mixing::Pulay_Mixing(this->mixing_ndim, this->mixing_beta);
            this->xdata.resize(this->mixing_ndim, 3, sizeof(double));
        }
        else if (method == "plain")
        {
            this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
            this->xdata.resize(1, 3, sizeof(double));
        }

        double* x_in = new double[3];
        this->x_out = new double[3];
        double* delta_x = new double[3];
        for (int i = 0; i < 3; ++i)
            x_in[i] = x_out[i] = 0.;

        auto screen = std::bind(&Mixing_Test::Kerker_mock, this, std::placeholders::_1);
        auto inner_dot = std::bind(&Mixing_Test::inner_dot_mock, this, std::placeholders::_1, std::placeholders::_2);
        
        double residual = 10.;

        while (niter < maxiter)
        {
            this->iteration(x_in, x_out);
            niter++;

            for (int i = 0; i < 3; ++i)
            {
                delta_x[i] = x_out[i] - x_in[i];
            }
            residual = this->inner_dot_mock(delta_x, delta_x);
            if (residual <= thr)
            {
                break;
            }

            this->mixing->push_data(this->xdata, x_in, x_out, screen, true);

            this->mixing->cal_coef(this->xdata, inner_dot);

            this->mixing->mix_data(this->xdata, x_in);
        }

        delete[] x_in;
        delete[] delta_x;
    }

    void iteration(double* x_in, double* x_out)
    {
        x_out[0] = (3. * x_in[1] - 2. * x_in[2] + 20.) / 8.;
        x_out[1] = (-4. * x_out[0] + 1. * x_in[2] + 33.) / 11.;
        x_out[2] = (-6. * x_out[0] - 3. * x_out[1] + 36.) / 12.;
    }

    void Kerker_mock(double *drho){}

    double inner_dot_mock(double* x1, double* x2)
    {
        double xnorm = 0.0;
        for (int ir = 0; ir < 3; ++ir)
        {
            xnorm += x1[ir] * x2[ir];
        }
        return xnorm;
    }
};


TEST_F(Mixing_Test, Broyden_Solve_LinearEq)
{
    solve("broyden");
    EXPECT_NEAR(x_out[0], 3., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1., DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 5);
}

TEST_F(Mixing_Test, Pulay_Solve_LinearEq)
{
    solve("pulay");
    EXPECT_NEAR(x_out[0], 3., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1., DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 5);
}

TEST_F(Mixing_Test, Plain_Solve_LinearEq)
{
    solve("plain");
    EXPECT_NEAR(x_out[0], 3.0000000013930039344, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 1.9999999779880637263, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1.0000000048064823233, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 11);
}