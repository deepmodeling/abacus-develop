#include "../symmetry_rotation_spin.h"

#include "source_base/constants.h"

#include <cmath>

#include "gtest/gtest.h"

using ModuleSymmetry::SpinRotation::Su2;
namespace SR = ModuleSymmetry::SpinRotation;

namespace
{
// Row-vector cartesian rotation gmatc for a column-vector rotation R_col:
//   gmatc = R_col^T (ABACUS row-vector convention).
// Column-vector rotation by angle `ang` about the z axis.
ModuleBase::Matrix3 col_rot_z(double ang)
{
    const double c = std::cos(ang), s = std::sin(ang);
    return ModuleBase::Matrix3(c, -s, 0, s, c, 0, 0, 0, 1);
}
ModuleBase::Matrix3 col_rot_x(double ang)
{
    const double c = std::cos(ang), s = std::sin(ang);
    return ModuleBase::Matrix3(1, 0, 0, 0, c, -s, 0, s, c);
}
ModuleBase::Matrix3 col_rot_y(double ang)
{
    const double c = std::cos(ang), s = std::sin(ang);
    return ModuleBase::Matrix3(c, 0, s, 0, 1, 0, -s, 0, c);
}

void expect_mat3_near(const ModuleBase::Matrix3& a, const ModuleBase::Matrix3& b, double tol = 1e-10)
{
    EXPECT_NEAR(a.e11, b.e11, tol);
    EXPECT_NEAR(a.e12, b.e12, tol);
    EXPECT_NEAR(a.e13, b.e13, tol);
    EXPECT_NEAR(a.e21, b.e21, tol);
    EXPECT_NEAR(a.e22, b.e22, tol);
    EXPECT_NEAR(a.e23, b.e23, tol);
    EXPECT_NEAR(a.e31, b.e31, tol);
    EXPECT_NEAR(a.e32, b.e32, tol);
    EXPECT_NEAR(a.e33, b.e33, tol);
}

void expect_su2_near(const Su2& a, const Su2& b, double tol = 1e-10)
{
    for (int i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(a[i].real(), b[i].real(), tol) << "elem " << i;
        EXPECT_NEAR(a[i].imag(), b[i].imag(), tol) << "elem " << i;
    }
}

bool is_unitary(const Su2& U, double tol = 1e-10)
{
    const Su2 prod = SR::mat2_mul(U, SR::dagger(U));
    return std::abs(prod[0] - 1.0) < tol && std::abs(prod[1]) < tol && std::abs(prod[2]) < tol
           && std::abs(prod[3] - 1.0) < tol;
}

std::complex<double> det2(const Su2& U)
{
    return U[0] * U[3] - U[1] * U[2];
}
} // namespace

// Identity operation -> U = I, W = I.
TEST(SymmetryRotationSpin, Identity)
{
    ModuleBase::Matrix3 g; // default ctor is identity
    const Su2 U = SR::so3_to_su2(g);
    expect_su2_near(U, Su2{1.0, 0.0, 0.0, 1.0});
    expect_mat3_near(SR::spin_so3(g), ModuleBase::Matrix3());
    expect_mat3_near(SR::pauli_rotation_matrix(U), ModuleBase::Matrix3());
}

// C4 about z (theta = pi/2). For a z-rotation U is diagonal diag(e^{-i th/2}, e^{i th/2}).
TEST(SymmetryRotationSpin, C4z)
{
    const double th = ModuleBase::PI / 2.0;
    const ModuleBase::Matrix3 gmatc = col_rot_z(th).Transpose(); // row-vector gmatc
    const Su2 U = SR::so3_to_su2(gmatc);
    const Su2 ref = {std::polar(1.0, -0.5 * th), 0.0, 0.0, std::polar(1.0, 0.5 * th)};
    expect_su2_near(U, ref);
    EXPECT_TRUE(is_unitary(U));
    EXPECT_NEAR(det2(U).real(), 1.0, 1e-12);
    EXPECT_NEAR(det2(U).imag(), 0.0, 1e-12);
}

// C2 about x: theta = pi special case. U(x, pi) = -i sigma_x = [[0,-i],[-i,0]].
TEST(SymmetryRotationSpin, C2x_ThetaPi)
{
    const ModuleBase::Matrix3 gmatc = col_rot_x(ModuleBase::PI).Transpose();
    const Su2 U = SR::so3_to_su2(gmatc);
    // up to double-group sign; fix sign by matching the (0,1) element direction
    Su2 ref = {std::complex<double>(0, 0), std::complex<double>(0, -1),
               std::complex<double>(0, -1), std::complex<double>(0, 0)};
    if ((U[1] + ref[1]).imag() == 0.0 && std::abs(U[1] - ref[1]) > 1e-6)
    {
        for (auto& z : ref)
        {
            z = -z;
        }
    }
    expect_su2_near(U, ref);
    EXPECT_TRUE(is_unitary(U));
    // W must be independent of the global U sign:
    expect_mat3_near(SR::pauli_rotation_matrix(U), SR::spin_so3(gmatc));
}

// Inversion: spin is a pseudovector -> proper part is identity -> U = I, W = I.
TEST(SymmetryRotationSpin, Inversion)
{
    const ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);
    const Su2 U = SR::so3_to_su2(inv);
    expect_su2_near(U, Su2{1.0, 0.0, 0.0, 1.0});
    expect_mat3_near(SR::spin_so3(inv), ModuleBase::Matrix3());
}

// Mirror plane z->-z (improper, det=-1): proper part is C2 about z.
TEST(SymmetryRotationSpin, MirrorZ)
{
    const ModuleBase::Matrix3 mz(1, 0, 0, 0, 1, 0, 0, 0, -1); // det = -1
    const Su2 U = SR::so3_to_su2(mz);
    EXPECT_TRUE(is_unitary(U));
    // W = diag(-1,-1,1): in-plane spin flips, out-of-plane preserved.
    expect_mat3_near(SR::pauli_rotation_matrix(U),
                     ModuleBase::Matrix3(-1, 0, 0, 0, -1, 0, 0, 0, 1));
    expect_mat3_near(SR::spin_so3(mz), ModuleBase::Matrix3(-1, 0, 0, 0, -1, 0, 0, 0, 1));
}

// Core consistency: for any operation, the SU(2) U built by so3_to_su2 induces the same
// Pauli (spin-vector) rotation as the closed-form W = R_proper^T. Sweep many angles/axes,
// proper and improper.
TEST(SymmetryRotationSpin, PauliConsistencySweep)
{
    std::vector<ModuleBase::Matrix3> cols;
    for (int k = 0; k <= 12; ++k)
    {
        const double a = ModuleBase::PI * k / 6.0;
        cols.push_back(col_rot_x(a));
        cols.push_back(col_rot_y(a));
        cols.push_back(col_rot_z(a));
    }
    // a few compound rotations
    cols.push_back(col_rot_z(0.7) * col_rot_y(1.3) * col_rot_x(2.1));
    cols.push_back(col_rot_x(2.5) * col_rot_z(1.1));

    for (const auto& Rcol : cols)
    {
        for (double det : {1.0, -1.0})
        {
            // build a row-vector gmatc, optionally improper (multiply by inversion)
            ModuleBase::Matrix3 gmatc = Rcol.Transpose();
            if (det < 0)
            {
                gmatc = gmatc * (-1.0);
            }
            const Su2 U = SR::so3_to_su2(gmatc);
            EXPECT_TRUE(is_unitary(U)) << "U not unitary";
            EXPECT_NEAR(det2(U).real(), 1.0, 1e-9);
            EXPECT_NEAR(det2(U).imag(), 0.0, 1e-9);
            // the two independent routes to W must agree
            expect_mat3_near(SR::pauli_rotation_matrix(U), SR::spin_so3(gmatc), 1e-9);
        }
    }
}

// Rotating a spin block U m U^dagger then by the inverse returns the original.
TEST(SymmetryRotationSpin, SpinBlockRoundTrip)
{
    const ModuleBase::Matrix3 g = col_rot_y(0.9).Transpose();
    const ModuleBase::Matrix3 ginv = g.Inverse();
    const Su2 U = SR::so3_to_su2(g);
    const Su2 Uinv = SR::so3_to_su2(ginv);
    const Su2 block = {std::complex<double>(0.3, 0.0), std::complex<double>(0.1, -0.2),
                       std::complex<double>(0.1, 0.2), std::complex<double>(-0.3, 0.0)};
    const Su2 rotated = SR::rotate_spin_block(block, U);
    const Su2 back = SR::rotate_spin_block(rotated, Uinv);
    expect_su2_near(back, block, 1e-9);
}

// Pauli-component rotation of a real-space density point matches W applied to (mx,my,mz).
TEST(SymmetryRotationSpin, RotatePauliComponents)
{
    const ModuleBase::Matrix3 g = col_rot_z(ModuleBase::PI / 3.0).Transpose();
    const double in[4] = {2.0, 1.0, 0.0, 0.5};
    double out[4];
    SR::rotate_pauli_components(g, in, out);
    EXPECT_NEAR(out[0], in[0], 1e-12); // charge untouched
    const ModuleBase::Matrix3 W = SR::spin_so3(g);
    EXPECT_NEAR(out[1], W.e11 * in[1] + W.e12 * in[2] + W.e13 * in[3], 1e-12);
    EXPECT_NEAR(out[2], W.e21 * in[1] + W.e22 * in[2] + W.e23 * in[3], 1e-12);
    EXPECT_NEAR(out[3], W.e31 * in[1] + W.e32 * in[2] + W.e33 * in[3], 1e-12);
    // magnitude of the spin vector is preserved by a proper rotation
    const double m2_in = in[1] * in[1] + in[2] * in[2] + in[3] * in[3];
    const double m2_out = out[1] * out[1] + out[2] * out[2] + out[3] * out[3];
    EXPECT_NEAR(m2_in, m2_out, 1e-10);
}
