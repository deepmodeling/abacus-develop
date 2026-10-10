#include "../exx_proj.h"

#include <gtest/gtest.h>
#include <vector>

namespace
{
const int naos = 5;
const int nocc = 2;
const int nvirt = 3;

std::vector<double> coefficients(const int columns, const double offset)
{
    std::vector<double> c(naos * columns);
    for (int col = 0; col < columns; ++col)
    {
        for (int row = 0; row < naos; ++row)
        {
            c[col * naos + row] = offset + 0.13 * row - 0.21 * col;
        }
    }
    return c;
}

// Independent reference: construct every rank-one probe density and contract
// it with a nonsymmetric AO matrix, as DMBand + cal_energy did.
std::vector<double> reference(const std::vector<double>& h,
                              const std::vector<double>& left,
                              const std::vector<double>& right,
                              const int nright,
                              const double factor)
{
    std::vector<double> result(nocc * nright, 0.0);
    for (int io = 0; io < nocc; ++io)
    {
        for (int iv = 0; iv < nright; ++iv)
        {
            for (int mu = 0; mu < naos; ++mu)
            {
                for (int nu = 0; nu < naos; ++nu)
                {
                    const double density = left[io * naos + mu] * right[iv * naos + nu];
                    result[io * nright + iv] += factor * density * h[nu * naos + mu];
                }
            }
        }
    }
    return result;
}

TEST(ExxProjection, AllFourModesPreserveNonsymmetricProbeContraction)
{
    std::vector<double> h(naos * naos);
    for (int col = 0; col < naos; ++col)
    {
        for (int row = 0; row < naos; ++row)
        {
            h[col * naos + row] = 0.17 * row + 0.32 * col + 0.09 * row * col;
        }
    }
    const auto co = coefficients(nocc, 0.3);
    const auto cv = coefficients(nvirt, -0.2);
    const auto cvx = coefficients(nocc, 0.7);
    const auto coxt = coefficients(nvirt, -0.5);
    const double factor = 0.5;
    const double minus_factor = -factor;
    for (int mode = 0; mode < 4; ++mode)
    {
        const int nright = mode == 1 ? nocc : nvirt;
        std::vector<double> result(nocc * nright, 1.25);
        std::vector<double> scratch(naos * nocc);
        std::vector<double> expected;
        if (mode == 0 || mode == 1)
        {
            const auto& right = mode == 1 ? co : cv;
            LR::project_exx(h.data(), co.data(), right.data(), naos, nocc, nright,
                            factor, scratch.data(), result.data());
            expected = reference(h, co, right, nright, factor);
        }
        else if (mode == 2)
        {
            LR::project_exx(h.data(), cvx.data(), cv.data(), naos, nocc, nvirt,
                            factor, scratch.data(), result.data());
            LR::project_exx(h.data(), co.data(), coxt.data(), naos, nocc, nvirt,
                            minus_factor, scratch.data(), result.data());
            expected = reference(h, cvx, cv, nvirt, factor);
            const auto occ = reference(h, co, coxt, nvirt, minus_factor);
            for (std::size_t i = 0; i < expected.size(); ++i) { expected[i] += occ[i]; }
        }
        else
        {
            LR::project_exx(h.data(), co.data(), coxt.data(), naos, nocc, nvirt,
                            factor, scratch.data(), result.data());
            expected = reference(h, co, coxt, nvirt, factor);
        }
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
            EXPECT_NEAR(result[i], 1.25 + expected[i], 1e-12) << "mode=" << mode << " element=" << i;
        }
    }
}

TEST(ExxProjection, EmptyProjectionDoesNotDereferenceInputs)
{
    double result = 3.0;
    LR::project_exx(nullptr, nullptr, nullptr, naos, 0, nvirt, 1.0, nullptr, &result);
    EXPECT_DOUBLE_EQ(result, 3.0);
}

TEST(ExxProjection, ComplexBlochPhasesAndAllFourProbeModes)
{
    typedef std::complex<double> Complex;
    const double factor = 0.5;
    const double minus_factor = -factor;
    const double angles[] = {0.0, 0.37, -0.61};
    std::vector<std::vector<Complex>> blocks(3, std::vector<Complex>(naos * naos));
    std::vector<Complex> co(naos * nocc);
    std::vector<Complex> cv(naos * nvirt);
    std::vector<Complex> cvx(naos * nocc);
    std::vector<Complex> coxt(naos * nvirt);
    for (int i = 0; i < naos * nvirt; ++i)
    {
        cv[i] = Complex(0.2 - 0.03 * i, 0.13 + 0.07 * i);
        coxt[i] = Complex(-0.1 + 0.02 * i, 0.04 - 0.05 * i);
    }
    for (int i = 0; i < naos * nocc; ++i)
    {
        co[i] = Complex(0.3 + 0.04 * i, -0.11 + 0.09 * i);
        cvx[i] = Complex(0.1 - 0.07 * i, 0.23 + 0.03 * i);
    }
    for (int r = 0; r < 3; ++r)
    {
        for (int i = 0; i < naos * naos; ++i)
        {
            blocks[r][i] = Complex(0.17 * i + 0.2 * r, 0.32 - 0.03 * i + 0.11 * r);
        }
    }
    // Multiple distinct k points and real-space cells. The reference forms
    // the original negative-phase probe and conjugates it in dotc(D,H).
    for (int ik = 0; ik < 3; ++ik)
    {
        std::vector<Complex> h(naos * naos, Complex(0));
        for (int r = 0; r < 3; ++r)
        {
            const Complex phase = std::exp(Complex(0, (ik + 1) * angles[r]));
            for (int i = 0; i < naos * naos; ++i) { h[i] += phase * blocks[r][i]; }
        }
        for (int mode = 0; mode < 4; ++mode)
        {
            const int nright = mode == 1 ? nocc : nvirt;
            const auto& left = mode == 2 ? cvx : co;
            const auto& right = mode == 1 ? co : (mode == 3 ? coxt : cv);
            const Complex initial(1.25, -0.3);
            std::vector<Complex> result(nocc * nright, initial);
            std::vector<Complex> scratch(naos * nocc);
            LR::project_exx(h.data(), left.data(), right.data(), naos, nocc, nright,
                            factor, scratch.data(), result.data());
            if (mode == 2)
            {
                LR::project_exx(h.data(), co.data(), coxt.data(), naos, nocc, nvirt,
                                minus_factor, scratch.data(), result.data());
            }
            for (int io = 0; io < nocc; ++io)
            {
                for (int iv = 0; iv < nright; ++iv)
                {
                    Complex expected = initial;
                    for (int r = 0; r < 3; ++r)
                    {
                        const Complex phase = std::exp(Complex(0, -(ik + 1) * angles[r]));
                        for (int mu = 0; mu < naos; ++mu)
                        {
                            for (int nu = 0; nu < naos; ++nu)
                            {
                                Complex density = phase * right[iv * naos + mu]
                                                * std::conj(left[io * naos + nu]);
                                if (mode == 1)
                                {
                                    density = phase * co[io * naos + mu] * std::conj(co[iv * naos + nu]);
                                }
                                if (mode == 2)
                                {
                                    density -= phase * coxt[iv * naos + mu] * std::conj(co[io * naos + nu]);
                                }
                                expected += factor * std::conj(density) * blocks[r][nu * naos + mu];
                            }
                        }
                    }
                    const int index = mode == 1 ? iv * nocc + io : io * nright + iv;
                    EXPECT_NEAR(std::abs(result[index] - expected), 0.0, 1e-11)
                        << "k=" << ik << " mode=" << mode << " io=" << io << " iv=" << iv;
                }
            }
        }
    }
}

TEST(ExxProjection, ComplexEmptyProjectionDoesNotDereferenceInputs)
{
    std::complex<double> result(3.0, 2.0);
    const std::complex<double>* input = nullptr;
    LR::project_exx(input, input, input, naos, nocc, 0, 1.0, nullptr, &result);
    EXPECT_EQ(result, std::complex<double>(3.0, 2.0));
}


TEST(ExxProjection, BenchmarkMappedIndicesPreserveProbeContraction)
{
    typedef std::complex<double> Complex;
    // Two atom types with 2 and 3 orbitals. The historical benchmark probe
    // uses those counts as coefficient indices for every orbital of a type.
    const int basis_index[] = {2, 2, 3, 3, 3};
    std::vector<double> original(naos * naos);
    std::vector<double> folded(naos * naos, 0.0);
    std::vector<Complex> original_complex(naos * naos);
    std::vector<Complex> folded_complex(naos * naos, Complex(0));
    for (int mu = 0; mu < naos; ++mu)
    {
        for (int nu = 0; nu < naos; ++nu)
        {
            const int index = nu * naos + mu;
            original[index] = 0.11 * mu - 0.23 * nu + 0.07 * mu * nu;
            original_complex[index] = Complex(original[index], 0.03 * mu + 0.13 * nu);
            const int mapped = basis_index[nu] * naos + basis_index[mu];
            folded[mapped] += original[index];
            folded_complex[mapped] += original_complex[index];
        }
    }
    const auto co = coefficients(nocc, 0.3);
    const auto cv = coefficients(nvirt, -0.2);
    std::vector<Complex> co_complex(co.begin(), co.end());
    std::vector<Complex> cv_complex(cv.begin(), cv.end());
    for (int i = 0; i < naos * nocc; ++i) { co_complex[i] += Complex(0, 0.05 * i); }
    for (int i = 0; i < naos * nvirt; ++i) { cv_complex[i] += Complex(0, -0.09 * i); }
    std::vector<double> result(nocc * nvirt, 0.0);
    std::vector<double> scratch(naos * nocc);
    std::vector<Complex> result_complex(nocc * nvirt, Complex(0));
    std::vector<Complex> scratch_complex(naos * nocc);
    const double factor = 0.5;
    LR::project_exx(folded.data(), co.data(), cv.data(), naos, nocc, nvirt,
                    factor, scratch.data(), result.data());
    LR::project_exx(folded_complex.data(), co_complex.data(), cv_complex.data(),
                    naos, nocc, nvirt, factor, scratch_complex.data(), result_complex.data());
    for (int io = 0; io < nocc; ++io)
    {
        for (int iv = 0; iv < nvirt; ++iv)
        {
            double expected = 0.0;
            Complex expected_complex(0);
            for (int mu = 0; mu < naos; ++mu)
            {
                for (int nu = 0; nu < naos; ++nu)
                {
                    const int row = basis_index[mu];
                    const int col = basis_index[nu];
                    expected += factor * co[io * naos + row] * cv[iv * naos + col]
                              * original[nu * naos + mu];
                    const Complex probe = cv_complex[iv * naos + row]
                                        * std::conj(co_complex[io * naos + col]);
                    expected_complex += factor * std::conj(probe) * original_complex[nu * naos + mu];
                }
            }
            const int index = io * nvirt + iv;
            EXPECT_NEAR(result[index], expected, 1e-12);
            EXPECT_NEAR(std::abs(result_complex[index] - expected_complex), 0.0, 1e-12);
        }
    }
}


TEST(ExxProjection, TransposedResponseMatchesSwappedCxcOProbe)
{
    std::vector<double> h(naos * naos);
    std::vector<double> transposed(naos * naos);
    for (int mu = 0; mu < naos; ++mu)
    {
        for (int nu = 0; nu < naos; ++nu)
        {
            const double value = 0.17 * mu - 0.32 * nu + 0.09 * mu * nu;
            h[nu * naos + mu] = value;
            transposed[mu * naos + nu] = value;
        }
    }
    const auto co = coefficients(nocc, 0.3);
    const auto coxt = coefficients(nvirt, -0.5);
    const double factor = 0.5;
    std::vector<double> result(nocc * nvirt, 0.0);
    std::vector<double> scratch(naos * nocc);
    LR::project_exx(transposed.data(), co.data(), coxt.data(), naos, nocc, nvirt,
                    factor, scratch.data(), result.data());
    const auto normal = reference(h, co, coxt, nvirt, factor);
    double difference = 0.0;
    for (int io = 0; io < nocc; ++io)
    {
        for (int iv = 0; iv < nvirt; ++iv)
        {
            double expected = 0.0;
            for (int mu = 0; mu < naos; ++mu)
            {
                for (int nu = 0; nu < naos; ++nu)
                {
                    expected += factor * coxt[iv * naos + mu] * co[io * naos + nu]
                              * h[nu * naos + mu];
                }
            }
            const int index = io * nvirt + iv;
            EXPECT_NEAR(result[index], expected, 1e-12);
            difference += std::abs(result[index] - normal[index]);
        }
    }
    EXPECT_GT(difference, 1e-6);
}

}
