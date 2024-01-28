#include "../inverse_matrix.h"
#include "gtest/gtest.h"
#include <iostream>

/************************************************
 *  unit test of inverse_matrix.h
 ***********************************************/

/**
 * - Tested Functions:
 *   - InverseMatrixComplex
 *     - use Inverse_Matrix_Complex to inverse a Hermite matrix
 *     - functions: init and using_zheev
 *
 *   - Inverse_Matrix_Real
 *     - computes the inverse of a dim*dim real matrix
 */

//a mock function of WARNING_QUIT, to avoid the uncorrected call by matrix.cpp at line 37.
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {exit(1);}
}

TEST(InverseMatrixComplexTest, InverseMatrixComplex)
{
    int dim = 10;
    ModuleBase::ComplexMatrix B(dim, dim);
    ModuleBase::ComplexMatrix C(dim, dim);
    ModuleBase::ComplexMatrix D(dim, dim);
    double a;
    double b;
    double c;
    // construct a Hermite matrix
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i <= j; i++)
        {
            if (i == j)
            {
                c = std::rand();
                B(i, j) = std::complex<double>(c, 0.0);
            }
            else
            {
                a = std::rand();
                b = std::rand();
                B(i, j) = std::complex<double>(a, b);
                B(j, i) = conj(B(i, j));
            }
        }
    }
    ModuleBase::Inverse_Matrix_Complex IMC;
    IMC.init(dim);
    IMC.using_zheev(B, C);
    D = B * C;
    for (int i = 0; i < dim; i++)
    {
        EXPECT_NEAR(D(i, i).real(), 1.0, 1e-14);
        EXPECT_NEAR(D(i, i).imag(), 0.0, 1e-14);
        // std::cout << D(i,i).real() << " " << D(i,i).imag()  << std::endl;
    }
}

TEST(InverseMatrixRealTest, InverseMatrixReal)
{
    int dim = 3;
    double in[9];
    double out[9];
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
            {
                in[i * dim + j] = 2.0;
            }
            else
            {
                in[i * dim + j] = 0.0;
            }
        }
    }
    ModuleBase::Inverse_Matrix_Real(dim, in, out);
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (i == j)
            {
                EXPECT_DOUBLE_EQ(in[i * dim + j], 2.0);
            }
            else
            {
                EXPECT_DOUBLE_EQ(in[i * dim + j], 0.0);
            }
        }
    }
    EXPECT_DOUBLE_EQ(out[0], 0.5);
    EXPECT_DOUBLE_EQ(out[1], 0.0);
    EXPECT_DOUBLE_EQ(out[2], 0.0);
    EXPECT_DOUBLE_EQ(out[3], 0.0);
    EXPECT_DOUBLE_EQ(out[4], 0.5);
    EXPECT_DOUBLE_EQ(out[5], 0.0);
    EXPECT_DOUBLE_EQ(out[6], 0.0);
    EXPECT_DOUBLE_EQ(out[7], 0.0);
    EXPECT_DOUBLE_EQ(out[8], 0.5);
}
