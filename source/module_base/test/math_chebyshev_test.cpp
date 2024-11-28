#include "../math_chebyshev.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of class Chebyshev
 ***********************************************/

/**
 * - Tested Functions:
 *   - calcoef_real
 *   - calcoef_complex
 *   - calcoef_pair
 *   - calfinalvec_real
 *   - calfinalvec_complex
 *   - tracepolyA
 *   - checkconverge
 *
 *
 */
class toolfunc
{
  public:
    double x7(double x)
    {
        return pow(x, 7);
    }
    double x6(double x)
    {
        return pow(x, 6);
    }
    double expr(double x)
    {
        return exp(x);
    }
    std::complex<double> expi(std::complex<double> x)
    {
        const std::complex<double> j(0.0, 1.0);
        return exp(j * x);
    }
    std::complex<double> expi2(std::complex<double> x)
    {
        const std::complex<double> j(0.0, 1.0);
        const double PI = 3.14159265358979323846;
        return exp(j * PI / 2.0 * x);
    }
    // Pauli matrix: [0,-i;i,0]
    int LDA = 2;
    double factor = 1;
    void sigma_y(std::complex<double>* spin_in, std::complex<double>* spin_out, const int m = 1)
    {
        const std::complex<double> j(0.0, 1.0);
        if (this->LDA < 2) {
            this->LDA = 2;
}
        for (int i = 0; i < m; ++i)
        {
            spin_out[LDA * i] = -factor * j * spin_in[LDA * i + 1];
            spin_out[LDA * i + 1] = factor * j * spin_in[LDA * i];
        }
    }
#ifdef __ENABLE_FLOAT_FFTW
    float x7(float x)
    {
        return pow(x, 7);
    }
    float x6(float x)
    {
        return pow(x, 6);
    }
    float expr(float x)
    {
        return exp(x);
    }
    std::complex<float> expi(std::complex<float> x)
    {
        const std::complex<float> j(0.0, 1.0);
        return exp(j * x);
    }
    std::complex<float> expi2(std::complex<float> x)
    {
        const std::complex<float> j(0.0, 1.0);
        const float PI = 3.14159265358979323846;
        return exp(j * PI / 2.0f * x);
    }
    // Pauli matrix: [0,-i;i,0]
    void sigma_y(std::complex<float>* spin_in, std::complex<float>* spin_out, const int m = 1)
    {
        const std::complex<float> j(0.0, 1.0);
        if (this->LDA < 2)
            this->LDA = 2;
        for (int i = 0; i < m; ++i)
        {
            spin_out[LDA * i] = -j * spin_in[LDA * i + 1];
            spin_out[LDA * i + 1] = j * spin_in[LDA * i];
        }
    }
#endif
};
class MathChebyshevTest : public testing::Test
{
  protected:
    ModuleBase::Chebyshev<double>* p_chetest;
    ModuleBase::Chebyshev<float>* p_fchetest;
    toolfunc fun;
};

TEST_F(MathChebyshevTest, calcoef_real)
{
    auto fun_x6 = [&](double x) { return fun.x6(x); };
    auto fun_x7 = [&](double x) { return fun.x7(x); };
    p_chetest = new ModuleBase::Chebyshev<double>(10);
    // x^6 = 1/32*( 10T_0 + 15T_2 + 6T_4 + T_6 )
    // x^7 = 1/64*( 35T_1 + 21T_3 + 7T_5 + T_7 )
    const double x6ref[10] = {10, 0, 15, 0, 6, 0, 1, 0, 0, 0};
    const double x7ref[10] = {0, 35, 0, 21, 0, 7, 0, 1, 0, 0};

    p_chetest->calcoef_real(fun_x6);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_chetest->coef_real[i] * 32.0, x6ref[i], 1.e-8);
    }
    p_chetest->calcoef_real(fun_x7);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_chetest->coef_real[i] * 64.0, x7ref[i], 1.e-8);
    }
    delete p_chetest;
}

TEST_F(MathChebyshevTest, calcoef_pair)
{
    auto fun_x6 = [&](double x) { return fun.x6(x); };
    auto fun_x7 = [&](double x) { return fun.x7(x); };
    p_chetest = new ModuleBase::Chebyshev<double>(10);
    // x^6 = 1/32*( 10T_0 + 15T_2 + 6T_4 + T_6 )
    // x^7 = 1/64*( 35T_1 + 21T_3 + 7T_5 + T_7 )
    const double x6ref[10] = {10, 0, 15, 0, 6, 0, 1, 0, 0, 0};
    const double x7ref[10] = {0, 35, 0, 21, 0, 7, 0, 1, 0, 0};
    p_chetest->calcoef_pair(fun_x6, fun_x7);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_chetest->coef_complex[i].real() * 32.0, x6ref[i], 1.e-8);
        EXPECT_NEAR(p_chetest->coef_complex[i].imag() * 64.0, x7ref[i], 1.e-8);
    }
    delete p_chetest;
}

TEST_F(MathChebyshevTest, calcoef_complex)
{
    auto fun_expi = [&](std::complex<double> x) { return fun.expi(x); };
    const int norder = 100;
    const double PI = 3.14159265358979323846;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    double* T = new double[norder];
    // check exp(i\pi/4) = \sum_n C_n[exp(ix)]T_n(\pi/4) = sqrt(2)/2*(1, i)
    p_chetest->calcoef_complex(fun_expi);
    p_chetest->getpolyval(PI / 4, T, norder);
    std::complex<double> sum(0, 0);
    for (int i = 0; i < norder; ++i)
    {
        sum += p_chetest->coef_complex[i] * T[i];
    }
    EXPECT_NEAR(sum.real(), sqrt(2) / 2, 1.e-8);
    EXPECT_NEAR(sum.imag(), sqrt(2) / 2, 1.e-8);
    delete[] T;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, calfinalvec_real)
{
    const int norder = 100;
    const double E = 2.718281828459046;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    //                 1  [ 1/e+e           -i(e-1/e) ]
    //  exp(\sigma_y)= -  [                           ], where \sigma_y = [0, -i; i, 0]
    //                 2  [ i(e-1/e)          1/e+e   ]
    std::complex<double>* v = new std::complex<double>[4];
    std::complex<double>* vout = new std::complex<double>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_expr = [&](double x) { return fun.expr(x); };
    auto fun_sigma_y
        = [&](std::complex<double>* in, std::complex<double>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    p_chetest->calcoef_real(fun_expr);
    p_chetest->calfinalvec_real(fun_sigma_y, v, vout, 2, 2, 2);
    EXPECT_NEAR(vout[0].real(), 0.5 * (E + 1 / E), 1.e-8);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[1].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[1].imag(), 0.5 * (E - 1 / E), 1.e-8);
    EXPECT_NEAR(vout[2].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[2].imag(), -0.5 * (E - 1 / E), 1.e-8);
    EXPECT_NEAR(vout[3].real(), 0.5 * (E + 1 / E), 1.e-8);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-8);

    delete[] v;
    delete[] vout;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, calfinalvec_complex)
{
    const int norder = 100;
    const double E = 2.718281828459046;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    //                        [ 0           1 ]
    //  exp(i pi/2*\sigma_y)= [               ], where \sigma_y = [0, -i; i, 0]
    //                        [ -1          0 ]
    std::complex<double>* v = new std::complex<double>[4];
    std::complex<double>* vout = new std::complex<double>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_sigma_y
        = [&](std::complex<double>* in, std::complex<double>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    auto fun_expi2 = [&](std::complex<double> x) { return fun.expi2(x); };
    p_chetest->calcoef_complex(fun_expi2);
    p_chetest->calfinalvec_complex(fun_sigma_y, v, vout, 2, 2, 2);
    EXPECT_NEAR(vout[0].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[1].real(), -1, 1.e-8);
    EXPECT_NEAR(vout[1].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[2].real(), 1, 1.e-8);
    EXPECT_NEAR(vout[2].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[3].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-8);

    delete[] v;
    delete[] vout;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, calpolyvec_complex)
{
    const int norder = 100;
    const double E = 2.718281828459046;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    //                        [ 0           1 ]
    //  exp(i pi/2*\sigma_y)= [               ], where \sigma_y = [0, -i; i, 0]
    //                        [ -1          0 ]
    std::complex<double>* v = new std::complex<double>[4];
    std::complex<double>* polyv = new std::complex<double>[4 * norder];
    std::complex<double>* vout = new std::complex<double>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]
    vout[0] = 0;
    vout[1] = 0;
    vout[2] = 0;
    vout[3] = 0;

    auto fun_sigma_y
        = [&](std::complex<double>* in, std::complex<double>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    auto fun_expi2 = [&](std::complex<double> x) { return fun.expi2(x); };
    p_chetest->calcoef_complex(fun_expi2);
    p_chetest->calpolyvec_complex(fun_sigma_y, v, polyv, 2, 2, 2);
    for (int i = 0; i < norder; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            vout[j] += polyv[i * 4 + j] * p_chetest->coef_complex[i];
        }
    }
    EXPECT_NEAR(vout[0].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[1].real(), -1, 1.e-8);
    EXPECT_NEAR(vout[1].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[2].real(), 1, 1.e-8);
    EXPECT_NEAR(vout[2].imag(), 0, 1.e-8);
    EXPECT_NEAR(vout[3].real(), 0, 1.e-8);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-8);

    delete[] v;
    delete[] vout;
    delete[] polyv;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, tracepolyA)
{
    const int norder = 100;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);

    // N == LDA
    std::complex<double>* v = new std::complex<double>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_sigma_y
        = [&](std::complex<double>* in, std::complex<double>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    p_chetest->tracepolyA(fun_sigma_y, v, 2, 2, 2);
    // Trace:  even function: 2 ; odd function 0.
    for (int i = 0; i < norder; ++i)
    {
        if (i % 2 == 0) {
            EXPECT_NEAR(p_chetest->polytrace[i], 2, 1.e-8);
        } else {
            EXPECT_NEAR(p_chetest->polytrace[i], 0, 1.e-8);
}
    }
    delete[] v;

    // N < LDA
    fun.LDA = 3;
    int LDA = fun.LDA;
    v = new std::complex<double>[2 * LDA];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 100.0;
    v[3] = 0.0;
    v[4] = 1.0;
    v[5] = 1.0; //[1 0; 0 1; 100 2]

    p_chetest->tracepolyA(fun_sigma_y, v, 2, LDA, 2);
    // Trace:  even function: 2 ; odd function 0.
    for (int i = 0; i < norder; ++i)
    {
        if (i % 2 == 0) {
            EXPECT_NEAR(p_chetest->polytrace[i], 2, 1.e-8);
        } else {
            EXPECT_NEAR(p_chetest->polytrace[i], 0, 1.e-8);
}
    }
    fun.LDA = 2;
    delete[] v;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, checkconverge)
{
    const int norder = 100;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    auto fun_sigma_y
        = [&](std::complex<double>* in, std::complex<double>* out, const int m = 1) { fun.sigma_y(in, out, m); };

    std::complex<double>* v = new std::complex<double>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]
    double tmin = -1.1;
    double tmax = 1.1;
    bool converge;
    converge = p_chetest->checkconverge(fun_sigma_y, v, 2, 2, tmax, tmin, 0.2);
    EXPECT_TRUE(converge);
    converge = p_chetest->checkconverge(fun_sigma_y, v + 2, 2, 2, tmax, tmin, 0.2);
    EXPECT_TRUE(converge);
    EXPECT_NEAR(tmin, -1.1, 1e-8);
    EXPECT_NEAR(tmax, 1.1, 1e-8);

    tmax = -1.1;
    converge = p_chetest->checkconverge(fun_sigma_y, v, 2, 2, tmax, tmin, 2.2);
    EXPECT_TRUE(converge);
    EXPECT_NEAR(tmin, -1.1, 1e-8);
    EXPECT_NEAR(tmax, 1.1, 1e-8);

    // not converge
    v[0] = std::complex<double>(0, 1), v[1] = 1;
    fun.factor = 1.5;
    tmin = -1.1, tmax = 1.1;
    converge = p_chetest->checkconverge(fun_sigma_y, v, 2, 2, tmax, tmin, 0.2);
    EXPECT_FALSE(converge);

    fun.factor = -1.5;
    tmin = -1.1, tmax = 1.1;
    converge = p_chetest->checkconverge(fun_sigma_y, v, 2, 2, tmax, tmin, 0.2);
    EXPECT_FALSE(converge);
    fun.factor = 1;

    delete[] v;
    delete p_chetest;
}

TEST_F(MathChebyshevTest, recurs)
{
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleBase::Chebyshev<double> noneche(0), ::testing::ExitedWithCode(1), "");
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("NOTICE"));

    int norder = 100;
    p_chetest = new ModuleBase::Chebyshev<double>(norder);
    int Tnp = p_chetest->recurs(1, 1, 0);
    EXPECT_EQ(Tnp, 2);
    delete p_chetest;
}

#ifdef __ENABLE_FLOAT_FFTW
TEST_F(MathChebyshevTest, calcoef_real_float)
{
    auto fun_x6f = [&](float x) { return fun.x6(x); };
    auto fun_x7f = [&](float x) { return fun.x7(x); };
    p_fchetest = new ModuleBase::Chebyshev<float>(10);
    // x^6 = 1/32*( 10T_0 + 15T_2 + 6T_4 + T_6 )
    // x^7 = 1/64*( 35T_1 + 21T_3 + 7T_5 + T_7 )
    const float x6ref[10] = {10, 0, 15, 0, 6, 0, 1, 0, 0, 0};
    const float x7ref[10] = {0, 35, 0, 21, 0, 7, 0, 1, 0, 0};
    p_fchetest->calcoef_real(fun_x6f);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_fchetest->coef_real[i] * 32.0, x6ref[i], 1.e-5);
    }
    p_fchetest->calcoef_real(fun_x7f);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_fchetest->coef_real[i] * 64.0, x7ref[i], 1.e-5);
    }
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, calcoef_pair_float)
{
    auto fun_x6f = [&](float x) { return fun.x6(x); };
    auto fun_x7f = [&](float x) { return fun.x7(x); };
    p_fchetest = new ModuleBase::Chebyshev<float>(10);
    // x^6 = 1/32*( 10T_0 + 15T_2 + 6T_4 + T_6 )
    // x^7 = 1/64*( 35T_1 + 21T_3 + 7T_5 + T_7 )
    const float x6ref[10] = {10, 0, 15, 0, 6, 0, 1, 0, 0, 0};
    const float x7ref[10] = {0, 35, 0, 21, 0, 7, 0, 1, 0, 0};
    p_fchetest->calcoef_pair(fun_x6f, fun_x7f);
    for (int i = 0; i < 10; ++i)
    {
        EXPECT_NEAR(p_fchetest->coef_complex[i].real() * 32.0, x6ref[i], 1.e-5);
        EXPECT_NEAR(p_fchetest->coef_complex[i].imag() * 64.0, x7ref[i], 1.e-5);
    }
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, calcoef_complex_float)
{
    auto fun_expif = [&](std::complex<float> x) { return fun.expi(x); };
    const int norder = 100;
    const float PI = 3.14159265358979323846;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);
    float* T = new float[norder];
    // check exp(i\pi/4) = \sum_n C_n[exp(ix)]T_n(\pi/4) = sqrt(2)/2*(1, i)
    p_fchetest->calcoef_complex(fun_expif);
    p_fchetest->getpolyval(PI / 4, T, norder);
    std::complex<float> sum(0, 0);
    for (int i = 0; i < norder; ++i)
    {
        sum += p_fchetest->coef_complex[i] * T[i];
    }
    EXPECT_NEAR(sum.real(), sqrt(2) / 2, 1.e-6);
    EXPECT_NEAR(sum.imag(), sqrt(2) / 2, 1.e-6);
    delete[] T;
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, calfinalvec_real_float)
{
    const int norder = 100;
    const float E = 2.718281828459046;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);
    //                 1  [ 1/e+e           -i(e-1/e) ]
    //  exp(\sigma_y)= -  [                           ], where \sigma_y = [0, -i; i, 0]
    //                 2  [ i(e-1/e)          1/e+e   ]
    std::complex<float>* v = new std::complex<float>[4];
    std::complex<float>* vout = new std::complex<float>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_exprf = [&](float x) { return fun.expr(x); };
    auto fun_sigma_yf
        = [&](std::complex<float>* in, std::complex<float>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    p_fchetest->calcoef_real(fun_exprf);
    p_fchetest->calfinalvec_real(fun_sigma_yf, v, vout, 2, 2, 2);
    EXPECT_NEAR(vout[0].real(), 0.5 * (E + 1 / E), 1.e-6);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[1].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[1].imag(), 0.5 * (E - 1 / E), 1.e-6);
    EXPECT_NEAR(vout[2].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[2].imag(), -0.5 * (E - 1 / E), 1.e-6);
    EXPECT_NEAR(vout[3].real(), 0.5 * (E + 1 / E), 1.e-6);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-6);

    delete[] v;
    delete[] vout;
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, calfinalvec_complex_float)
{
    const int norder = 100;
    const float E = 2.718281828459046;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);
    //                        [ 0           1 ]
    //  exp(i pi/2*\sigma_y)= [               ], where \sigma_y = [0, -i; i, 0]
    //                        [ -1          0 ]
    std::complex<float>* v = new std::complex<float>[4];
    std::complex<float>* vout = new std::complex<float>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_sigma_yf
        = [&](std::complex<float>* in, std::complex<float>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    auto fun_expi2f = [&](std::complex<float> x) { return fun.expi2(x); };
    p_fchetest->calcoef_complex(fun_expi2f);
    p_fchetest->calfinalvec_complex(fun_sigma_yf, v, vout, 2, 2, 2);
    EXPECT_NEAR(vout[0].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[1].real(), -1, 1.e-6);
    EXPECT_NEAR(vout[1].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[2].real(), 1, 1.e-6);
    EXPECT_NEAR(vout[2].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[3].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-6);

    delete[] v;
    delete[] vout;
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, calpolyvec_float)
{
    const int norder = 100;
    const float E = 2.718281828459046;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);
    //                        [ 0           1 ]
    //  exp(i pi/2*\sigma_y)= [               ], where \sigma_y = [0, -i; i, 0]
    //                        [ -1          0 ]
    std::complex<float>* v = new std::complex<float>[4];
    std::complex<float>* polyv = new std::complex<float>[4 * norder];
    std::complex<float>* vout = new std::complex<float>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]
    vout[0] = 0;
    vout[1] = 0;
    vout[2] = 0;
    vout[3] = 0;

    auto fun_sigma_yf
        = [&](std::complex<float>* in, std::complex<float>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    auto fun_expi2f = [&](std::complex<float> x) { return fun.expi2(x); };
    p_fchetest->calcoef_complex(fun_expi2f);
    p_fchetest->calpolyvec_complex(fun_sigma_yf, v, polyv, 2, 2, 2);
    for (int i = 0; i < norder; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            vout[j] += polyv[i * 4 + j] * p_fchetest->coef_complex[i];
        }
    }
    EXPECT_NEAR(vout[0].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[0].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[1].real(), -1, 1.e-6);
    EXPECT_NEAR(vout[1].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[2].real(), 1, 1.e-6);
    EXPECT_NEAR(vout[2].imag(), 0, 1.e-6);
    EXPECT_NEAR(vout[3].real(), 0, 1.e-6);
    EXPECT_NEAR(vout[3].imag(), 0, 1.e-6);

    delete[] v;
    delete[] vout;
    delete[] polyv;
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, tracepolyA_float)
{
    const int norder = 100;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);

    std::complex<float>* v = new std::complex<float>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]

    auto fun_sigma_yf
        = [&](std::complex<float>* in, std::complex<float>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    p_fchetest->tracepolyA(fun_sigma_yf, v, 2, 2, 2);
    // Trace:  even function: 2 ; odd function 0.
    for (int i = 0; i < norder; ++i)
    {
        if (i % 2 == 0)
            EXPECT_NEAR(p_fchetest->polytrace[i], 2, 1.e-6);
        else
            EXPECT_NEAR(p_fchetest->polytrace[i], 0, 1.e-6);
    }
    delete[] v;

    // N < LDA
    fun.LDA = 3;
    int LDA = fun.LDA;
    v = new std::complex<float>[2 * LDA];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 100.0;
    v[3] = 0.0;
    v[4] = 1.0;
    v[5] = 1.0; //[1 0; 0 1; 100 2]

    p_fchetest->tracepolyA(fun_sigma_yf, v, 2, LDA, 2);
    // Trace:  even function: 2 ; odd function 0.
    for (int i = 0; i < norder; ++i)
    {
        if (i % 2 == 0)
            EXPECT_NEAR(p_fchetest->polytrace[i], 2, 1.e-6);
        else
            EXPECT_NEAR(p_fchetest->polytrace[i], 0, 1.e-6);
    }
    fun.LDA = 2;
    delete[] v;
    delete p_fchetest;
}

TEST_F(MathChebyshevTest, checkconverge_float)
{
    const int norder = 100;
    p_fchetest = new ModuleBase::Chebyshev<float>(norder);

    std::complex<float>* v = new std::complex<float>[4];
    v[0] = 1.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 1.0; //[1 0; 0 1]
    float tmin = -1.1;
    float tmax = 1.1;
    bool converge;

    auto fun_sigma_yf
        = [&](std::complex<float>* in, std::complex<float>* out, const int m = 1) { fun.sigma_y(in, out, m); };
    converge = p_fchetest->checkconverge(fun_sigma_yf, v, 2, 2, tmax, tmin, 0.2);
    EXPECT_TRUE(converge);
    converge = p_fchetest->checkconverge(fun_sigma_yf, v + 2, 2, 2, tmax, tmin, 0.2);
    EXPECT_TRUE(converge);
    EXPECT_NEAR(tmin, -1.1, 1e-6);
    EXPECT_NEAR(tmax, 1.1, 1e-6);

    delete[] v;
    delete p_fchetest;
}
#endif