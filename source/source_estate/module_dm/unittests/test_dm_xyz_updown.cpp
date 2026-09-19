#include <cmath>
#include <complex>

#include "gtest/gtest.h"

#include "source_estate/module_dm/density_matrix.h"

// The runtime stores D_ab = conj(c_a) c_b = conj(P_ab).
// These tests recover the physical Pauli traces from that stored block,
// including nonzero y moments, both template types, and nonzero offsets.
// PR #7832 corrects the old test's confusion between stored D and physical P.

namespace
{
constexpr double TOL = 1e-12;
// contiguous 2x2 spin block at icol = 0 (as in cal_DMR for a single orbital pair)
const int step_trace[4] = {0, 1, 2, 3};

void expect_dvec_near(const double* got, const double r0, const double rx, const double ry, const double rz)
{
    EXPECT_NEAR(got[0], r0, TOL);
    EXPECT_NEAR(got[1], rx, TOL);
    EXPECT_NEAR(got[2], ry, TOL);
    EXPECT_NEAR(got[3], rz, TOL);
}

void expect_cvec_near(const std::complex<double>* got, const std::complex<double>* ref)
{
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(got[i].real(), ref[i].real(), TOL) << "elem " << i;
        EXPECT_NEAR(got[i].imag(), ref[i].imag(), TOL) << "elem " << i;
    }
}
} // namespace

// Pure spinor psi = (1+i, 2-0.5i): D = conj(psi psi^dagger) is Hermitian, so the
// complex-template output must be purely real and identical to the double
// template. Reference values from m_i = Tr(sigma_i rho):
//   m_x = 2 Re(psi_up^* psi_dn) = 3, m_y = 2 Im(psi_up^* psi_dn) = -5,
//   m_z = |psi_up|^2 - |psi_dn|^2 = -2.25, charge = |psi|^2 = 6.25.
TEST(FuncXyzToUpdown, HermitianSpinor)
{
    const std::complex<double> up(1.0, 1.0);
    const std::complex<double> dn(2.0, -0.5);
    const std::complex<double> tmp[4] = {std::norm(up), std::conj(up) * dn, up * std::conj(dn), std::norm(dn)};
    // Stored D_ud = conj(1+i)*(2-0.5i) = 1.5 - 2.5i.
    EXPECT_NEAR(tmp[1].real(), 1.5, TOL);
    EXPECT_NEAR(tmp[1].imag(), -2.5, TOL);

    double out_d[4];
    std::complex<double> out_c[4];
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);

    expect_dvec_near(out_d, 6.25, 3.0, -5.0, -2.25);
    const std::complex<double> ref[4] = {6.25, 3.0, -5.0, -2.25};
    expect_cvec_near(out_c, ref);
}

// Axis magnetizations from the actual conjugate-first storage. Both templates agree.
TEST(FuncXyzToUpdown, AxisMagnetizations)
{
    // spin +x: psi = (1,1)/sqrt2
    {
        const std::complex<double> up = std::complex<double>(1.0, 0.0) / std::sqrt(2.0);
        const std::complex<double> dn = std::complex<double>(1.0, 0.0) / std::sqrt(2.0);
        const std::complex<double> tmp[4] = {std::norm(up), std::conj(up) * dn, up * std::conj(dn), std::norm(dn)};
        double out_d[4];
        std::complex<double> out_c[4];
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);
        expect_dvec_near(out_d, 1.0, 1.0, 0.0, 0.0);
        const std::complex<double> ref[4] = {1.0, 1.0, 0.0, 0.0};
        expect_cvec_near(out_c, ref);
    }
    // spin +y: psi = (1,i)/sqrt2  ->  rho_y = +1
    {
        const std::complex<double> up = std::complex<double>(1.0, 0.0) / std::sqrt(2.0);
        const std::complex<double> dn = std::complex<double>(0.0, 1.0) / std::sqrt(2.0);
        const std::complex<double> tmp[4] = {std::norm(up), std::conj(up) * dn, up * std::conj(dn), std::norm(dn)};
        double out_d[4];
        std::complex<double> out_c[4];
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);
        expect_dvec_near(out_d, 1.0, 0.0, 1.0, 0.0);
        const std::complex<double> ref[4] = {1.0, 0.0, 1.0, 0.0};
        expect_cvec_near(out_c, ref);
    }
    // spin -y: psi = (1,-i)/sqrt2 ->  rho_y = -1
    {
        const std::complex<double> up = std::complex<double>(1.0, 0.0) / std::sqrt(2.0);
        const std::complex<double> dn = std::complex<double>(0.0, -1.0) / std::sqrt(2.0);
        const std::complex<double> tmp[4] = {std::norm(up), std::conj(up) * dn, up * std::conj(dn), std::norm(dn)};
        double out_d[4];
        std::complex<double> out_c[4];
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);
        expect_dvec_near(out_d, 1.0, 0.0, -1.0, 0.0);
        const std::complex<double> ref[4] = {1.0, 0.0, -1.0, 0.0};
        expect_cvec_near(out_c, ref);
    }
}

// Generic (not necessarily Hermitian) 4-tuple: pins the exact formulas of both
// templates, including that the double template equals the real part of the
// complex template componentwise.
TEST(FuncXyzToUpdown, GenericTuple)
{
    const std::complex<double> tmp[4] = {std::complex<double>(1.0, 2.0),
                                         std::complex<double>(0.5, 0.3),
                                         std::complex<double>(-0.7, 0.1),
                                         std::complex<double>(0.3, -0.4)};
    double out_d[4];
    std::complex<double> out_c[4];
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);

    // double: rho_0 = 1+0.3, rho_x = 0.5-0.7, rho_y = 0.3-0.1, rho_z = 1-0.3
    expect_dvec_near(out_d, 1.3, -0.2, 0.2, 0.7);
    // complex: rho_y = -i*(tmp1 - tmp2) = -i*(1.2+0.2i) = 0.2-1.2i
    const std::complex<double> ref[4] = {std::complex<double>(1.3, 1.6),
                                         std::complex<double>(-0.2, 0.4),
                                         std::complex<double>(0.2, -1.2),
                                         std::complex<double>(0.7, 2.4)};
    expect_cvec_near(out_c, ref);
    // cross-consistency: the double template is the real part of the complex one
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(out_d[i], out_c[i].real(), TOL) << "elem " << i;
    }
}

// Pauli -> spinor -> Pauli round trip: build the spinor block from a given
// (rho_0, rho_x, rho_y, rho_z), stored as D_ud=(rho_x+i*rho_y)/2, and check that
// func_xyz_to_updown recovers the original components.
TEST(FuncXyzToUpdown, PauliRoundTrip)
{
    const double r0 = 2.0, rx = -1.5, ry = 0.75, rz = 0.25;
    const std::complex<double> tmp[4] = {std::complex<double>(0.5 * (r0 + rz), 0.0),
                                         std::complex<double>(0.5 * rx, 0.5 * ry),
                                         std::complex<double>(0.5 * rx, -0.5 * ry),
                                         std::complex<double>(0.5 * (r0 - rz), 0.0)};
    double out_d[4];
    std::complex<double> out_c[4];
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out_d);
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<std::complex<double>>(tmp, 0, step_trace, out_c);
    expect_dvec_near(out_d, r0, rx, ry, rz);
    const std::complex<double> ref[4] = {r0, rx, ry, rz};
    expect_cvec_near(out_c, ref);
}

// icol/step_trace offset handling: with a nonzero base column the outputs land
// at target[icol + step_trace[k]] (the same indexing used by cal_DMR).
TEST(FuncXyzToUpdown, IndexOffset)
{
    const std::complex<double> tmp[4] = {std::complex<double>(1.0, 0.0),
                                         std::complex<double>(0.0, 1.0),
                                         std::complex<double>(0.0, -1.0),
                                         std::complex<double>(2.0, 0.0)};
    const int col_size = 8; // arbitrary column stride of the atom-pair block
    const int step[4] = {0, 1, col_size, col_size + 1};
    double out[2 * col_size + 2] = {0.0};
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, col_size / 2, step, out);
    // rho_0 = 3, rho_x = 0, rho_y = 1 - (-1) = 2, rho_z = -1 at icol = 4
    EXPECT_NEAR(out[col_size / 2 + 0], 3.0, TOL);
    EXPECT_NEAR(out[col_size / 2 + 1], 0.0, TOL);
    EXPECT_NEAR(out[col_size / 2 + col_size], 2.0, TOL);
    EXPECT_NEAR(out[col_size / 2 + col_size + 1], -1.0, TOL);
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
