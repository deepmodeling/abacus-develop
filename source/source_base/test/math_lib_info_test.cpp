#define GATHER_INFO
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/lapack_connector.h"
#undef GATHER_INFO

#include "gtest/gtest.h"
#include <complex>

TEST(MathLibInfoTest, DelegatesBlasOperations)
{
    const char no_transpose = 'N';
    const int one = 1;
    const int two = 2;
    const std::complex<double> alpha(2.0, 0.0);
    const std::complex<double> beta(0.0, 0.0);
    const std::complex<double> a(3.0, 0.0);
    const std::complex<double> b(4.0, 0.0);
    std::complex<double> c(0.0, 0.0);

    zgemm_(&no_transpose, &no_transpose, &one, &one, &one, &alpha, &a, &one, &b, &one, &beta, &c, &one);
    EXPECT_EQ(c, std::complex<double>(24.0, 0.0));

    const std::complex<double> x[2] = {{1.0, 1.0}, {2.0, -1.0}};
    std::complex<double> y[2] = {{0.0, 0.0}, {1.0, 1.0}};
    zaxpy_(&two, &alpha, x, &one, y, &one);
    EXPECT_EQ(y[0], std::complex<double>(2.0, 2.0));
    EXPECT_EQ(y[1], std::complex<double>(5.0, -1.0));
}

TEST(MathLibInfoTest, DelegatesGeneralizedEigenproblem)
{
    const int problem_type = 1;
    const char eigenvectors = 'V';
    const char all_eigenvalues = 'A';
    const char upper_triangle = 'U';
    const int one = 1;
    const double lower_bound = 0.0;
    const double upper_bound = 0.0;
    const double absolute_tolerance = 0.0;
    const int workspace_size = 4;
    std::complex<double> a[1] = {{4.0, 0.0}};
    std::complex<double> b[1] = {{2.0, 0.0}};
    int eigenvalue_count = 0;
    double eigenvalue[1] = {0.0};
    std::complex<double> eigenvector[1] = {{0.0, 0.0}};
    std::complex<double> workspace[workspace_size];
    double real_workspace[7] = {0.0};
    int integer_workspace[5] = {0};
    int failed_indices[1] = {0};
    int info = -1;

    zhegvx_(&problem_type,
            &eigenvectors,
            &all_eigenvalues,
            &upper_triangle,
            &one,
            a,
            &one,
            b,
            &one,
            &lower_bound,
            &upper_bound,
            &one,
            &one,
            &absolute_tolerance,
            &eigenvalue_count,
            eigenvalue,
            eigenvector,
            &one,
            workspace,
            &workspace_size,
            real_workspace,
            integer_workspace,
            failed_indices,
            &info);

    EXPECT_EQ(info, 0);
    EXPECT_EQ(eigenvalue_count, 1);
    EXPECT_NEAR(eigenvalue[0], 2.0, 1.0e-14);
}
