#include "gtest/gtest.h"
#include "../broyden_mixing.h"

#define DOUBLETHRESHOLD 1e-8

class Broyden_Mixing_Test : public testing::Test
{
protected:
    const double mixing_beta = 0.8;
    const int mixing_ndim = 7;
    Base_Mixing::Mixing_Data xdata;
    Base_Mixing::Broyden_Mixing* broyden = nullptr;

    double thr = 1e-5;
    int niter = 0;
    int maxiter = 50;
    double* x_out = nullptr;

    void SetUp()
    {
        this->broyden = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
    }
    void TearDown()
    {
        delete this->broyden;
        delete[] this->x_out;
    }

    void solve()
    {
        double* x_in = new double[3];
        this->x_out = new double[3];
        double* delta_x = new double[3];
        for (int i = 0; i < 3; ++i)
            x_in[i] = x_out[i] = 0.;

        auto screen = std::bind(&Broyden_Mixing_Test::Kerker_mock, this, std::placeholders::_1);
        auto inner_dot = std::bind(&Broyden_Mixing_Test::inner_dot_mock, this, std::placeholders::_1, std::placeholders::_2);
        this->xdata.resize(this->mixing_ndim, 3, sizeof(double));
        
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

            this->broyden->push_data(this->xdata, x_in, x_out, screen, true);

            this->broyden->cal_coef(this->xdata, inner_dot);

            this->broyden->mix_data(this->xdata, x_in);
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


TEST_F(Broyden_Mixing_Test, Solve_LinearEq)
{
    solve();
    EXPECT_NEAR(x_out[0], 3., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2., DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1., DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 5);
}
