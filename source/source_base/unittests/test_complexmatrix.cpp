#include <complex>
#include <stdexcept>

#include <gtest/gtest.h>

#include "source_base/complexmatrix.h"

TEST(ComplexMatrix, ElementAccessUsesSizeTIndex)
{
    ModuleBase::ComplexMatrix m(3, 4, true);
    m(2, 3) = std::complex<double>(1.0, -2.0);
    EXPECT_DOUBLE_EQ(std::real(m(2, 3)), 1.0);
    EXPECT_DOUBLE_EQ(std::imag(m(2, 3)), -2.0);
}

TEST(ComplexMatrix, SizeIsElementCount)
{
    ModuleBase::ComplexMatrix m(3, 4);
    EXPECT_EQ(m.size, 12u);
}

TEST(ComplexMatrix, NegativeConstructorThrows)
{
    EXPECT_THROW(ModuleBase::ComplexMatrix m(-1, 4), std::length_error);
}

TEST(ComplexMatrix, CreateNegativeThrows)
{
    ModuleBase::ComplexMatrix m;
    EXPECT_THROW(m.create(-1, 4), std::length_error);
}

TEST(ComplexMatrix, EmptyMatrixHasNoStorage)
{
    ModuleBase::ComplexMatrix m(0, 4);
    EXPECT_EQ(m.nr, 0);
    EXPECT_EQ(m.nc, 4);
    EXPECT_EQ(m.size, 0u);
    EXPECT_EQ(m.c, nullptr);
}

TEST(ComplexMatrix, IdentityMatrix)
{
    ModuleBase::ComplexMatrix m(2, 2, false);
    m.set_as_identity_matrix();
    EXPECT_DOUBLE_EQ(std::real(m(0, 0)), 1.0);
    EXPECT_DOUBLE_EQ(std::real(m(1, 1)), 1.0);
    EXPECT_DOUBLE_EQ(std::real(m(0, 1)), 0.0);
    EXPECT_DOUBLE_EQ(std::real(m(1, 0)), 0.0);
}

TEST(ComplexMatrix, FromRealMatrix)
{
    ModuleBase::matrix r(3, 4);
    for (int i = 0; i < 12; ++i)
    {
        r.c[i] = static_cast<double>(i);
    }

    ModuleBase::ComplexMatrix m(r);
    EXPECT_EQ(m.size, 12u);
    EXPECT_DOUBLE_EQ(std::real(m(2, 3)), 11.0);
    EXPECT_DOUBLE_EQ(std::imag(m(2, 3)), 0.0);
}
