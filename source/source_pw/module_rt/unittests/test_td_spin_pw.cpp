#include "source_pw/module_rt/unittests/td_test_utils.h"

#include <cmath>
#include <limits>

namespace td_test
{
TEST_F(Observables, SOCTorqueMatchesCommutator)
{
    const ModuleBase::Vector3<double> a(0.13, -0.07, 0.09);
    const std::string mag_partition = "none";
    pw::TDSpinPW<double, Device> observer;
    observer.initialize(cell_, basis_, basis_, kv_, sf_, pp_, mag_partition, "", 0.0, true);
    const pw::TDSpinResultPW actual = observer.calculate(*psi_, state_, a);
    const Vector beta = projectors(a);
    const int n = psi_->get_nbasis();
    double signal = 0;
    for (int axis = 0; axis < 3; ++axis)
    {
        double total = 0;
        for (int atom = 0; atom < cell_.nat; ++atom)
        {
            double expected = 0;
            for (int band = 0; band < psi_->get_nbands(); ++band)
            {
                const Complex* data = psi_->get_pointer() + band * n;
                const Vector x(data, data + n);
                const Vector sx = pauli(x, axis);
                const Vector hx = nonlocal(beta, x, atom);
                const Vector hsx = nonlocal(beta, sx, atom);
                const Vector shx = pauli(hx, axis);
                const Complex commutator = inner(x, hsx) - inner(x, shx);
                expected += state_.wg(0, band) * std::real(Complex(0, 0.5) * commutator);
            }
            EXPECT_NEAR(actual.atom_soc[3 * atom + axis], expected, 1e-10);
            total += expected;
            signal += std::abs(expected);
        }
        EXPECT_NEAR(actual.soc[axis], total, 1e-10);
    }
    EXPECT_GT(signal, 1e-8);
}

TEST_F(Observables, DensityMomentsAndPeriodicVoronoi)
{
    const std::vector<int> labels = pw::td_voronoi_pw(cell_, basis_);
    std::vector<double> atom_integrals(cell_.nat, 0.0);
    const double dv = cell_.omega / basis_.nxyz;
    for (int ir = 0; ir < basis_.nrxx; ++ir)
    {
        const int ix = ir / (basis_.ny * basis_.nplane);
        const int iy = ir / basis_.nplane % basis_.ny;
        const int iz = ir % basis_.nplane + basis_.startz_current;
        const ModuleBase::Vector3<double> point(double(ix) / basis_.nx, double(iy) / basis_.ny, double(iz) / basis_.nz);
        double best = std::numeric_limits<double>::infinity();
        int nearest = -1;
        // Independent exhaustive periodic-image reference for this moderately skew cell.
        for (int atom = 0; atom < cell_.nat; ++atom)
        {
            for (int i = -3; i <= 3; ++i)
            {
                for (int j = -3; j <= 3; ++j)
                {
                    for (int k = -3; k <= 3; ++k)
                    {
                        const ModuleBase::Vector3<double> image(i, j, k);
                        const ModuleBase::Vector3<double> delta = (point - cell_.atoms[0].taud[atom] + image) * cell_.latvec;
                        const double distance = delta.norm2();
                        if (distance < best - 1e-12)
                        {
                            best = distance;
                            nearest = atom;
                        }
                    }
                }
            }
        }
        ASSERT_EQ(labels[ir], nearest);
        const double f = 1.0 + 0.01 * ir;
        atom_integrals[nearest] += f * dv;
        charge_.rho[0][ir] = 2 * f;
        charge_.rho[1][ir] = 0.2 * f;
        charge_.rho[2][ir] = -0.3 * f;
        charge_.rho[3][ir] = 0.5 * f;
        // Potential storage is in Ry; the reference spin field is (0.4, 0.1, -0.2) Ha.
        state_.pot->get_eff_v()(1, ir) = 0.8;
        state_.pot->get_eff_v()(2, ir) = 0.2;
        state_.pot->get_eff_v()(3, ir) = -0.4;
    }
    const std::string mag_partition = "voronoi";
    pw::TDSpinPW<double, Device> observer;
    observer.initialize(cell_, basis_, basis_, kv_, sf_, pp_, mag_partition, "", 0.0, true);
    const ModuleBase::Vector3<double> zero(0, 0, 0);
    const pw::TDSpinResultPW actual = observer.calculate(*psi_, state_, zero);
    const double integral = cell_.omega * (1.0 + 0.005 * (basis_.nxyz - 1));
    EXPECT_NEAR(actual.electrons, 2 * integral, 1e-10);
    EXPECT_NEAR(actual.mag_abs, std::sqrt(0.38) * integral, 1e-10);
    const double m[3] = {0.2, -0.3, 0.5};
    const double torque[3] = {-0.01, -0.24, -0.14};
    for (int axis = 0; axis < 3; ++axis)
    {
        EXPECT_NEAR(actual.mag[axis], m[axis] * integral, 1e-10);
        EXPECT_NEAR(actual.xc[axis], torque[axis] * integral, 1e-10);
        double sum = 0;
        for (int atom = 0; atom < cell_.nat; ++atom)
        {
            EXPECT_NEAR(actual.atom_mag[3 * atom + axis], m[axis] * atom_integrals[atom], 1e-10);
            sum += actual.atom_mag[3 * atom + axis];
        }
        EXPECT_NEAR(sum, actual.mag[axis], 1e-10);
    }
}
} // namespace td_test
