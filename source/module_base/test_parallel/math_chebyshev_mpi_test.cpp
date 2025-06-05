#include "../math_chebyshev.h"
#include "mpi.h"
#include "module_base/parallel_comm.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
/************************************************
 *  unit test of class Chebyshev MPI part
 ***********************************************/

 /**
  * - Tested Functions:
  * - checkconverge
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
    int dsize = 0;
    int my_rank = 0;
    void SetUp() override
    {
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
        int color = (world_rank < world_size / 2) ? 0 : 1;
        int key = world_rank;
    
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &POOL_WORLD);
    
        int pool_rank, pool_size;
        MPI_Comm_rank(POOL_WORLD, &pool_rank);
        MPI_Comm_size(POOL_WORLD, &pool_size);
    }
    void TearDown() override
    {
    }
};

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

#ifdef __ENABLE_FLOAT_FFTW
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

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
