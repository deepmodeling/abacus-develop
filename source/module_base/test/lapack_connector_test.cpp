#include "../lapack_connector.h"

#include <cmath>
#include <complex>
#include <vector>

#include "gtest/gtest.h"

/************************************************
 *  unit test of lapack_connector.h
 ***********************************************/

/**
 * - Tested Functions:
 *  - zhegv_
 *   - use zhegv_ to compute the eigenvalues and eigenvectors of
 *   - a complex Hermitian-definite generalized eigenproblem
 */

class LapackConnectorTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        // Initialize matrices A and B and the eigenvalue vector
        // (Use appropriate values for your test case)
        A = {
            std::complex<double>(2.0, 0.0),
            std::complex<double>(1.0, -1.0),
            std::complex<double>(1.0, 1.0),
            std::complex<double>(3.0, 0.0),
        };

        // Create a random square matrix C with complex elements
        std::vector<std::complex<double>> C = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0},
            {7.0, 8.0}
        };

        // Compute the conjugate transpose of C
        std::vector<std::complex<double>> C_conj_transpose = {
            {C[0].real(), -C[0].imag()},
            {C[1].real(), -C[1].imag()},
            {C[2].real(), -C[2].imag()},
            {C[3].real(), -C[3].imag()}
        };

        // Compute the product of C_conj_transpose and C to obtain B
        B = {{C_conj_transpose[0] * C[0] + C_conj_transpose[1] * C[1]},
             {C_conj_transpose[0] * C[2] + C_conj_transpose[1] * C[3]},
             {C_conj_transpose[2] * C[0] + C_conj_transpose[3] * C[1]},
             {C_conj_transpose[2] * C[2] + C_conj_transpose[3] * C[3]}};

        n = sqrt(A.size());
        lda = n;
        ldb = n;
        w.resize(n);

        // Set up the parameters for zhegv_
        itype = 1;
        jobz = 'V';
        uplo = 'U';
        lwork = -1;
        info = 0;

        // Ensure that B is positive definite
    }

    int itype;
    char jobz;
    char uplo;
    int n;
    int lda;
    int ldb;
    int lwork;
    int info;
    // matrices A and B are column-major
    std::vector<std::complex<double>> A;
    std::vector<std::complex<double>> B;
    std::vector<double> w;
};

// Test the zhegv_ function
TEST_F(LapackConnectorTest, ZHEGV)
{
    // First, query the optimal size of the work array
    std::complex<double> work_query;
    double rwork_query;
    zhegv_(&itype,
           &jobz,
           &uplo,
           &n,
           A.data(),
           &lda,
           B.data(),
           &ldb,
           w.data(),
           &work_query,
           &lwork,
           &rwork_query,
           &info);
    lwork = static_cast<int>(work_query.real());
    std::vector<std::complex<double>> work(lwork);
    // std::vector<double> rwork(static_cast<int>(rwork_query));
    // the above line is not working as rwork_query will return -nan
    // std::vector<double> rwork(7 * lwork);
    std::vector<double> rwork(7 * n);

    // Now, call zhegv_ with the optimal work array size
    zhegv_(&itype,
           &jobz,
           &uplo,
           &n,
           A.data(),
           &lda,
           B.data(),
           &ldb,
           w.data(),
           work.data(),
           &lwork,
           rwork.data(),
           &info);

    // Check that the function completed successfully
    ASSERT_EQ(info, 0);

    // Check the computed eigenvalues and eigenvectors
    // (Use appropriate values for your test case)
    std::vector<double> expected_eigenvalues = {0.014371905048252809, 1.0871905949517402};
    std::vector<std::complex<double>> expected_eigenvectors = {
        {0.00029066041795582461, -0.042636598658647745},
        {0.07557994526773984,    0.0                  },
        {-0.81903769393029213,   -0.083945171943878405},
        {0.33387897788468901,    0.0                  }
    };

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(w[i], expected_eigenvalues[i], 1e-8);
        for (size_t j = 0; j < n; ++j)
        {
            EXPECT_NEAR(A[i * n + j].real(), expected_eigenvectors[i * n + j].real(), 1e-8);
            EXPECT_NEAR(A[i * n + j].imag(), expected_eigenvectors[i * n + j].imag(), 1e-8);
        }
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}