#include "source_pw/module_pwdft/kernels/stress_op.h"

#include <gtest/gtest.h>
#include <complex>
#include <vector>

namespace
{
// The DeltaSpin branch of cal_stress_nl_op accumulates, for one strain component,
//     stress += -1/2 * d(lambda . M)/d(eps)
// where dbecp is d(becp)/d(eps) for that component and M = b^dagger * sigma * b is
// the physical spinor moment. The 1/2 is the symmetric-strain convention shared by
// all three branches of this kernel: the DFT+U and collinear branches likewise use
// fac = d_wg without the explicit 2.0 that the force kernel carries. The fixtures
// below therefore pin the relative signs of the three Pauli channels against
// moments built independently of the kernel, and do not re-derive the overall
// strain normalization.
template <typename Real>
struct DeltaSpinStressFixture
{
    int nh[2] = {2, 1};
    int na[2] = {1, 1};
    int isk[2] = {0, 1};
    const int occupied = 2;
    const int nkb = 3;
    const int wg_nc = 4;
    const int ik = 1;
    std::vector<Real> weights = {0.9, 0.8, 0.7, 0.6, 0.3, 0.17, 0.0, 0.0};
    std::vector<Real> lambda;
    std::vector<std::complex<Real>> b;
    std::vector<std::complex<Real>> db;

    void prepare(const int npol, const int axis)
    {
        lambda.assign(6, 0);
        lambda[axis] = 0.7;
        lambda[3 + axis] = -0.4;
        b.resize(nkb * npol * occupied);
        db.resize(b.size());
        for (std::size_t i = 0; i < b.size(); ++i)
        {
            b[i] = std::complex<Real>(0.13 + 0.017 * i, -0.21 + 0.009 * i);
        }
        for (std::size_t i = 0; i < db.size(); ++i)
        {
            db[i] = std::complex<Real>(-0.07 + 0.003 * i, 0.11 - 0.005 * i);
        }
    }

    // Mirrors the kernel's projector to atom mapping (type 0: two projectors on the
    // first atom, type 1: one projector on the second), so only the spin algebra
    // differs between the oracle and the kernel.
    int atom_of(const int proj) const
    {
        int block = 0;
        int atom = 0;
        for (int it = 0; it < 2; ++it)
        {
            const int width = nh[it] * na[it];
            if (proj < block + width)
            {
                return atom + (proj - block) / nh[it];
            }
            block += width;
            atom += na[it];
        }
        return atom;
    }

    // Independent oracle: physical Pauli moments of b displaced along db.
    double energy(const int npol, const double displacement) const
    {
        double result = 0;
        const int sign = isk[ik] == 0 ? 1 : -1;
        for (int ib = 0; ib < occupied; ++ib)
        {
            const int row = npol == 2 ? ib * 2 * nkb : ib * nkb;
            const double weight = weights[ik * wg_nc + ib];
            for (int proj = 0; proj < nkb; ++proj)
            {
                const int index = row + proj;
                const int iat = atom_of(proj);
                const std::complex<double> up = std::complex<double>(b[index])
                    + displacement * std::complex<double>(db[index]);
                double value = sign * lambda[iat * 3 + 2] * std::norm(up);
                if (npol == 2)
                {
                    const std::complex<double> down = std::complex<double>(b[index + nkb])
                        + displacement * std::complex<double>(db[index + nkb]);
                    const std::complex<double> cross = std::conj(up) * down;
                    value = lambda[iat * 3] * 2 * cross.real()
                        + lambda[iat * 3 + 1] * 2 * cross.imag()
                        + lambda[iat * 3 + 2] * (std::norm(up) - std::norm(down));
                }
                result += weight * value;
            }
        }
        return result;
    }

    void check(const std::vector<Real>& stress, const int npol, const int axis, const double tolerance) const
    {
        const double h = 1e-4;
        const double ep = energy(npol, h);
        const double em = energy(npol, -h);
        const double expected = -0.5 * (ep - em) / (2 * h);
        // Start from a nonzero stress to also verify additive accumulation.
        EXPECT_NEAR(stress[0] - Real(0.125), expected, tolerance) << "npol=" << npol << " axis=" << axis;
    }

    void cpu(const int npol, std::vector<Real>& stress) const
    {
        hamilt::cal_stress_nl_op<Real, base_device::DEVICE_CPU>()(
            nullptr, nkb, occupied, 2, wg_nc, ik, npol, nh, na,
            weights.data(), lambda.data(), isk, b.data(), db.data(), stress.data());
    }
};

template <typename Real>
void check_cpu(const double tolerance)
{
    DeltaSpinStressFixture<Real> fixture;
    for (int axis = 0; axis < 3; ++axis)
    {
        fixture.prepare(2, axis);
        std::vector<Real> stress(1, 0.125);
        fixture.cpu(2, stress);
        fixture.check(stress, 2, axis, tolerance);
    }
    fixture.prepare(1, 2);
    for (int spin = 0; spin < 2; ++spin)
    {
        fixture.isk[1] = spin;
        std::vector<Real> stress(1, 0.125);
        fixture.cpu(1, stress);
        fixture.check(stress, 1, 2, tolerance);
    }
}

TEST(DeltaSpinStress, CpuDoubleEnergyDerivative)
{
    check_cpu<double>(1e-10);
}
} // namespace
