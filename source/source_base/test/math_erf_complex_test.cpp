#include "source_base/math_erf_complex.h"

#include "gtest/gtest.h"
#include <cmath>
#include <complex>

namespace
{

void expect_complex_near(const std::complex<double>& actual,
                         const std::complex<double>& expected,
                         const double tolerance)
{
    EXPECT_NEAR(actual.real(), expected.real(), tolerance);
    EXPECT_NEAR(actual.imag(), expected.imag(), tolerance);
}

} // namespace

TEST(ErrorFuncTest, RealArgumentsAgreeWithStandardLibrary)
{
    ModuleBase::ErrorFunc error_function;
    const double x = 0.75;

    expect_complex_near(ModuleBase::ErrorFunc::erf(std::complex<double>(x, 0.0)),
                        std::complex<double>(std::erf(x), 0.0),
                        1.0e-14);
    expect_complex_near(ModuleBase::ErrorFunc::erfc(std::complex<double>(x, 0.0)),
                        std::complex<double>(std::erfc(x), 0.0),
                        1.0e-14);
    EXPECT_NEAR(ModuleBase::ErrorFunc::erfcx(x), std::exp(x * x) * std::erfc(x), 1.0e-14);
}

TEST(ErrorFuncTest, ComplexFunctionsSatisfyDefiningIdentities)
{
    const std::complex<double> z(0.75, -0.5);
    const std::complex<double> imaginary_unit(0.0, 1.0);
    const std::complex<double> erf_z = ModuleBase::ErrorFunc::erf(z);
    const std::complex<double> erfc_z = ModuleBase::ErrorFunc::erfc(z);
    const std::complex<double> erfcx_z = ModuleBase::ErrorFunc::erfcx(z);
    const std::complex<double> erfi_z = ModuleBase::ErrorFunc::erfi(z);
    const std::complex<double> scaled_w_z = ModuleBase::ErrorFunc::scaled_w(z, 1.0e-13);

    expect_complex_near(erf_z + erfc_z, std::complex<double>(1.0, 0.0), 1.0e-13);
    expect_complex_near(erfcx_z, std::exp(z * z) * erfc_z, 1.0e-13);
    expect_complex_near(erfi_z, -imaginary_unit * ModuleBase::ErrorFunc::erf(imaginary_unit * z), 1.0e-13);
    expect_complex_near(scaled_w_z, std::exp(-z * z) * ModuleBase::ErrorFunc::erfc(-imaginary_unit * z), 1.0e-13);
}

TEST(ErrorFuncTest, RealFaddeevaImaginaryPartAndErfiAreConsistent)
{
    const double x = 1.25;
    const std::complex<double> scaled_w = ModuleBase::ErrorFunc::scaled_w(std::complex<double>(x, 0.0), 0.0);

    EXPECT_NEAR(scaled_w.real(), std::exp(-x * x), 1.0e-14);
    EXPECT_NEAR(scaled_w.imag(), ModuleBase::ErrorFunc::scaled_w_im(x), 1.0e-14);
    EXPECT_NEAR(ModuleBase::ErrorFunc::erfi(x), std::exp(x * x) * scaled_w.imag(), 1.0e-13);
}

TEST(ErrorFuncTest, SymmetryHoldsAcrossAlgorithmRegions)
{
    const std::complex<double> values[] = {
        std::complex<double>(1.0e-8, 2.0e-8),
        std::complex<double>(-1.0e-4, 1.0),
        std::complex<double>(3.0, 4.0),
        std::complex<double>(-8.0, 0.25),
        std::complex<double>(1.0, -12.0),
    };

    for (const std::complex<double>& z: values)
    {
        const std::complex<double> erf_z = ModuleBase::ErrorFunc::erf(z, 1.0e-12);
        const std::complex<double> erf_conjugate = ModuleBase::ErrorFunc::erf(std::conj(z), 1.0e-12);
        expect_complex_near(erf_conjugate, std::conj(erf_z), 1.0e-11);
    }
}
