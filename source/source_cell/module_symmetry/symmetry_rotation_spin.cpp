#include "symmetry_rotation_spin.h"

#include "source_base/constants.h"

#include <cmath>

namespace ModuleSymmetry
{
namespace SpinRotation
{
using cd = std::complex<double>;

namespace
{
// Proper part of an orthogonal cartesian operation: spin only sees the proper rotation.
// For det(gmatc) = +1 it is gmatc itself; for det = -1 (improper) it is -gmatc.
ModuleBase::Matrix3 proper_part(const ModuleBase::Matrix3& gmatc)
{
    return gmatc.Det() < 0.0 ? gmatc * (-1.0) : gmatc;
}

double clamp_cos(double c)
{
    if (c > 1.0)
    {
        return 1.0;
    }
    if (c < -1.0)
    {
        return -1.0;
    }
    return c;
}

// U(n, theta) = cos(theta/2) I - i sin(theta/2) (n . sigma)
Su2 su2_from_axis_angle(double nx, double ny, double nz, double theta)
{
    const double c = std::cos(0.5 * theta);
    const double s = std::sin(0.5 * theta);
    // -i s (n.sigma) =
    //   [ -i s nz          -i s nx - s ny ]
    //   [ -i s nx + s ny    i s nz        ]
    Su2 U;
    U[0] = cd(c, -s * nz);            // (uu)
    U[1] = cd(-s * ny, -s * nx);      // (ud)
    U[2] = cd(s * ny, -s * nx);       // (du)
    U[3] = cd(c, s * nz);             // (dd)
    return U;
}
} // namespace

Su2 so3_to_su2(const ModuleBase::Matrix3& gmatc, const double eps)
{
    const ModuleBase::Matrix3 R = proper_part(gmatc);

    const double trace = R.e11 + R.e22 + R.e33;
    const double cos_theta = clamp_cos(0.5 * (trace - 1.0));
    const double theta = std::acos(cos_theta);
    const double sin_theta = std::sin(theta);

    // theta ~ 0 : identity rotation, U = I
    if (theta < eps)
    {
        return Su2{cd(1.0, 0.0), cd(0.0, 0.0), cd(0.0, 0.0), cd(1.0, 0.0)};
    }

    // theta ~ pi : axis-angle formula is singular (sin theta -> 0). 
    // Extract the axis from the symmetric part: 
    // for theta = pi,  R = 2 n n^T - I, so n_i n_j = (R_ij + delta_ij)/2.
    if (std::abs(theta - ModuleBase::PI) < eps || std::abs(sin_theta) < eps)
    {
        double nn[3] = {0.5 * (R.e11 + 1.0), 0.5 * (R.e22 + 1.0), 0.5 * (R.e33 + 1.0)};
        for (int i = 0; i < 3; ++i)
        {
            nn[i] = nn[i] > 0.0 ? nn[i] : 0.0; // guard tiny negatives
        }
        // pick the largest diagonal as the reference component (sign fixed to +)
        int imax = 0;
        if (nn[1] > nn[imax])
        {
            imax = 1;
        }
        if (nn[2] > nn[imax])
        {
            imax = 2;
        }
        double n[3] = {0.0, 0.0, 0.0};
        n[imax] = std::sqrt(nn[imax]);
        // off-diagonal of (R + I)/2 gives n_i n_j -> fix relative signs from the reference
        const double off[3][3] = {{0.0, 0.5 * (R.e12 + R.e21), 0.5 * (R.e13 + R.e31)},
                                  {0.5 * (R.e21 + R.e12), 0.0, 0.5 * (R.e23 + R.e32)},
                                  {0.5 * (R.e31 + R.e13), 0.5 * (R.e32 + R.e23), 0.0}};
        for (int j = 0; j < 3; ++j)
        {
            if (j == imax)
            {
                continue;
            }
            n[j] = off[imax][j] / n[imax];
        }
        // normalize for safety
        const double norm = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        if (norm > 0.0)
        {
            n[0] /= norm;
            n[1] /= norm;
            n[2] /= norm;
        }
        return su2_from_axis_angle(n[0], n[1], n[2], ModuleBase::PI);
    }

    // general case: axis from the antisymmetric part, expressed in ROW-vector elements.
    // n_x = (R_yz - R_zy)/(2 sin),  n_y = (R_zx - R_xz)/(2 sin),  n_z = (R_xy - R_yx)/(2 sin)
    const double inv = 1.0 / (2.0 * sin_theta);
    const double nx = (R.e23 - R.e32) * inv;
    const double ny = (R.e31 - R.e13) * inv;
    const double nz = (R.e12 - R.e21) * inv;
    return su2_from_axis_angle(nx, ny, nz, theta);
}

ModuleBase::Matrix3 spin_so3(const ModuleBase::Matrix3& gmatc)
{
    // rho'^i = sum_j W_ij rho^j with W = R_proper^T (= column-vector rotation R_col).
    return proper_part(gmatc).Transpose();
}

ModuleBase::Matrix3 pauli_rotation_matrix(const Su2& U)
{
    // sigma matrices (row-major 2x2)
    static const Su2 sx = {cd(0, 0), cd(1, 0), cd(1, 0), cd(0, 0)};
    static const Su2 sy = {cd(0, 0), cd(0, -1), cd(0, 1), cd(0, 0)};
    static const Su2 sz = {cd(1, 0), cd(0, 0), cd(0, 0), cd(-1, 0)};
    const Su2 sig[3] = {sx, sy, sz};
    const Su2 Ud = dagger(U);

    double w[3][3];
    for (int j = 0; j < 3; ++j)
    {
        const Su2 rot = mat2_mul(mat2_mul(U, sig[j]), Ud); // U sigma_j U^dagger
        for (int i = 0; i < 3; ++i)
        {
            // W_ij = (1/2) Tr(sigma_i * rot)
            const Su2& si = sig[i];
            const cd tr = si[0] * rot[0] + si[1] * rot[2] + si[2] * rot[1] + si[3] * rot[3];
            w[i][j] = 0.5 * tr.real();
        }
    }
    return ModuleBase::Matrix3(w[0][0], w[0][1], w[0][2],
                               w[1][0], w[1][1], w[1][2],
                               w[2][0], w[2][1], w[2][2]);
}

Su2 dagger(const Su2& U)
{
    return Su2{std::conj(U[0]), std::conj(U[2]), std::conj(U[1]), std::conj(U[3])};
}

Su2 mat2_mul(const Su2& A, const Su2& B)
{
    return Su2{A[0] * B[0] + A[1] * B[2], A[0] * B[1] + A[1] * B[3],
               A[2] * B[0] + A[3] * B[2], A[2] * B[1] + A[3] * B[3]};
}

Su2 rotate_spin_block(const Su2& block, const Su2& U)
{
    return mat2_mul(mat2_mul(U, block), dagger(U));
}

void rotate_pauli_components(const ModuleBase::Matrix3& gmatc, const double in[4], double out[4])
{
    const ModuleBase::Matrix3 W = spin_so3(gmatc);
    out[0] = in[0]; // charge component is a scalar, untouched by spin rotation
    out[1] = W.e11 * in[1] + W.e12 * in[2] + W.e13 * in[3];
    out[2] = W.e21 * in[1] + W.e22 * in[2] + W.e23 * in[3];
    out[3] = W.e31 * in[1] + W.e32 * in[2] + W.e33 * in[3];
}
} // namespace SpinRotation
} // namespace ModuleSymmetry
