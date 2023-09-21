#include "../broyden_mixing.h"
#include "../plain_mixing.h"
#include "../pulay_mixing.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <omp.h>

#define DOUBLETHRESHOLD 1e-8
double ext_inner_product_mock(double* x1, double* x2)
{
}
class Mixing_Test : public testing::Test
{
  protected:
    Mixing_Test()
    {
    }
    ~Mixing_Test()
    {
        delete this->mixing;
    }
    const double mixing_beta = 0.6;
    const int mixing_ndim = 3;
    Base_Mixing::Mixing_Data xdata;
    Base_Mixing::Mixing* mixing = nullptr;
    double thr = 1e-8;
    int niter = 0;
    int maxiter = 10;
    std::vector<double> xd_ref = {0.0, 0.0, 0.0};
    std::vector<std::complex<double>> xc_ref = {
        {0.0, 1.0},
        {1.0, 0.0},
        0.0
    };
    void init_method(std::string method)
    {
        if (method == "broyden")
        {
            this->mixing = new Base_Mixing::Broyden_Mixing(this->mixing_ndim, this->mixing_beta);
        }
        else if (method == "pulay")
        {
            this->mixing = new Base_Mixing::Pulay_Mixing(this->mixing_ndim, this->mixing_beta);
        }
        else if (method == "plain")
        {
            this->mixing = new Base_Mixing::Plain_Mixing(this->mixing_beta);
        }
    }

    void clear()
    {
        delete this->mixing;
        this->mixing = nullptr;
    }

    /**
     * @brief sover linear equation:
     *        [ 8 -3  2 ][x1]   [20]       [3]
     *        [ 4 11 -1 ][x2] = [33]   x = [2]
     *        [ 6  3 12 ][x3]   [36]       [1]
     *
     *         [x1]   [ 3/8  -2/8   20/8 ][x1]
     *         [x2] = [-4/11  1/11 -33/11][x2]
     *         [x3]   [-6/12 -3/12  36/12][x3]
     */
    template <typename FPTYPE>
    void solve_linear_eq(FPTYPE* x_in, FPTYPE* x_out)
    {
        this->mixing->init_mixing_data(xdata, 3, sizeof(FPTYPE));
        std::vector<FPTYPE> delta_x(3);

        auto screen = std::bind(&Mixing_Test::Kerker_mock<FPTYPE>, this, std::placeholders::_1);
        auto inner_product
            = std::bind(static_cast<double (Mixing_Test::*)(FPTYPE*, FPTYPE*)>(&Mixing_Test::inner_product_mock),
                        this,
                        std::placeholders::_1,
                        std::placeholders::_2);

        double residual = 10.;
        this->niter = 0;
        while (niter < maxiter)
        {
            x_out[0] = (3. * x_in[1] - 2. * x_in[2] + 20.) / 8.;
            x_out[1] = (-4. * x_out[0] + 1. * x_in[2] + 33.) / 11.;
            x_out[2] = (-6. * x_out[0] - 3. * x_out[1] + 36.) / 12.;

            niter++;

            for (int i = 0; i < 3; ++i)
            {
                delta_x[i] = x_out[i] - x_in[i];
            }
            residual = this->inner_product_mock(delta_x.data(), delta_x.data());
            if (residual <= thr)
            {
                break;
            }

            this->mixing->push_data(this->xdata, x_in, x_out, screen, true);

            this->mixing->cal_coef(this->xdata, inner_product);

            this->mixing->mix_data(this->xdata, x_in);
        }
    }

    template <typename FPTYPE>
    void Kerker_mock(FPTYPE* drho)
    {
    }

    double inner_product_mock(double* x1, double* x2)
    {
        double xnorm = 0.0;
        for (int ir = 0; ir < 3; ++ir)
        {
            xnorm += x1[ir] * x2[ir];
        }
        return xnorm;
    }
    double inner_product_mock(std::complex<double>* x1, std::complex<double>* x2)
    {
        double xnorm = 0.0;
        for (int ir = 0; ir < 3; ++ir)
        {
            xnorm += x1[ir].real() * x2[ir].real() + x1[ir].imag() * x2[ir].imag();
        }
        return xnorm;
    }
};

TEST_F(Mixing_Test, Broyden_Solve_LinearEq)
{
    omp_set_num_threads(1);
    init_method("broyden");
    std::vector<double> x_in = xd_ref;
    std::vector<double> x_out(3);
    solve_linear_eq<double>(x_in.data(), x_out.data());
    EXPECT_NEAR(x_out[0], 2.9999999999999996, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2.0000000000000004, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1.0000000000000000, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 5);

    this->mixing->reset();
    xdata.reset();

    std::vector<std::complex<double>> xc_in = xc_ref;
    std::vector<std::complex<double>> xc_out(3);
    solve_linear_eq<std::complex<double>>(xc_in.data(), xc_out.data());
    EXPECT_NEAR(xc_out[0].real(), 3.0000000000000009, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[1].real(), 1.9999999999999998, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[2].real(), 0.99999999999999944, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 5);
    std::string output;
    Base_Mixing::Mixing_Data testdata;
    this->mixing->init_mixing_data(testdata, 3, sizeof(double));

    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->mixing->push_data(testdata, x_in.data(), x_out.data(), nullptr, true),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(
        output,
        testing::HasSubstr("One Broyden_Mixing object can only bind one Mixing_Data object to calculate coefficients"));

    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->mixing->cal_coef(testdata, ext_inner_product_mock),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(
        output,
        testing::HasSubstr("One Broyden_Mixing object can only bind one Mixing_Data object to calculate coefficients"));

    clear();
}

TEST_F(Mixing_Test, Pulay_Solve_LinearEq)
{
    omp_set_num_threads(1);
    init_method("pulay");
    std::vector<double> x_in = xd_ref;
    std::vector<double> x_out(3);
    solve_linear_eq<double>(x_in.data(), x_out.data());
    EXPECT_NEAR(x_out[0], 2.9999959638248037, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2.0000002552633349, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1.0000019542717642, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 6);

    this->mixing->reset();
    xdata.reset();

    std::vector<std::complex<double>> xc_in = xc_ref;
    std::vector<std::complex<double>> xc_out(3);
    solve_linear_eq<std::complex<double>>(xc_in.data(), xc_out.data());
    EXPECT_NEAR(xc_out[0].real(), 3.0000063220482565, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[1].real(), 1.9999939191147462, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[2].real(), 0.99999835919718549, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 6);

    std::string output;
    Base_Mixing::Mixing_Data testdata;
    this->mixing->init_mixing_data(testdata, 3, sizeof(double));

    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->mixing->push_data(testdata, x_in.data(), x_out.data(), nullptr, true),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(
        output,
        testing::HasSubstr("One Pulay_Mixing object can only bind one Mixing_Data object to calculate coefficients"));

    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->mixing->cal_coef(testdata, ext_inner_product_mock),
                ::testing::ExitedWithCode(0),
                "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(
        output,
        testing::HasSubstr("One Pulay_Mixing object can only bind one Mixing_Data object to calculate coefficients"));

    clear();
}

TEST_F(Mixing_Test, Plain_Solve_LinearEq)
{
    omp_set_num_threads(1);
    init_method("plain");
    std::vector<double> x_in = xd_ref;
    std::vector<double> x_out(3);
    solve_linear_eq<double>(x_in.data(), x_out.data());
    EXPECT_NEAR(x_out[0], 2.9999613068687698, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[1], 2.0000472873362103, DOUBLETHRESHOLD);
    EXPECT_NEAR(x_out[2], 1.0000075247315625, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 10);

    this->mixing->reset();
    xdata.reset();

    std::vector<std::complex<double>> xc_in = xc_ref;
    std::vector<std::complex<double>> xc_out(3);
    solve_linear_eq<std::complex<double>>(xc_in.data(), xc_out.data());
    EXPECT_NEAR(xc_out[0].real(), 2.9999418982632711, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[1].real(), 2.0000317031363761, DOUBLETHRESHOLD);
    EXPECT_NEAR(xc_out[2].real(), 1.0000211250842703, DOUBLETHRESHOLD);
    ASSERT_EQ(niter, 10);

    this->mixing->reset();

    clear();
}