#include "gtest/gtest.h"
#include "source_estate/module_dm/density_matrix.h"

#include <complex>
#include <cmath>

/************************************************************************
 * Regression test for the nspin=4 (non-collinear/SOC) magnetization
 * round-trip through the density-matrix pipeline.
 *
 * Physical invariant (must hold regardless of internal sign conventions):
 *   the magnetization <sigma> of the occupied one-electron state that is
 *   encoded in the density matrix must be recovered, with the CORRECT SIGN
 *   in ALL THREE cartesian components, by func_xyz_to_updown().
 *
 * Why this test exists (see ref/MSG/2026-08-#7664-...):
 *   ABACUS builds the k-space DM as  DM_{ab} = sum_n w_n conj(c_{n,a}) c_{n,b}
 *   (cal_dm_psi.cpp: the conj() is applied to the FIRST index a). Hence the
 *   stored DM block is the complex conjugate of the physical 1-RDM P:
 *       DM_{up,dn} = conj(c_up) c_dn = conj(P_{up,dn}).
 *   Since m_x, m_z read Re() (conjugation-invariant) but m_y reads Im(),
 *   ONLY m_y is sensitive to this conjugation. func_xyz_to_updown() must be
 *   consistent with that stored convention. PR #7664 set the m_y extraction
 *   to the "bare" textbook formula (valid for P, not for conj(P)), which
 *   flips m_y for in-plane moments and quenches non-collinear order
 *   (e.g. Mn3Sn 120-degree AFM). This test pins m_y down.
 *
 * The helper build_DM_block_as_cal_dm_psi() MUST mirror cal_dm_psi.cpp. If
 * that convention is ever changed (e.g. the "upstream" fix that makes the DM
 * hold the physical P), update the helper in the SAME commit so this test
 * keeps asserting the physical invariant.
 ************************************************************************/

namespace
{
using cd = std::complex<double>;

// spinor of the occupied state with <sigma> = mhat (the +1 eigenstate of mhat.sigma)
void spinor_from_direction(const double mhat[3], cd c[2])
{
    // |+n> = (cos(th/2), sin(th/2) e^{i ph});  n=(sin th cos ph, sin th sin ph, cos th)
    const double th = std::acos(std::max(-1.0, std::min(1.0, mhat[2])));
    const double ph = std::atan2(mhat[1], mhat[0]);
    c[0] = cd(std::cos(0.5 * th), 0.0);
    c[1] = std::sin(0.5 * th) * cd(std::cos(ph), std::sin(ph));
}

// Build the 4 spinor-block DM elements EXACTLY as cal_dm_psi.cpp stores them:
//   DM_{a,b} = sum_occ w * conj(c_a) * c_b     (conj on the first index)
// layout tmp = {uu, ud, du, dd}
void build_DM_block_as_cal_dm_psi(const cd c[2], double w, cd tmp[4])
{
    tmp[0] = w * std::conj(c[0]) * c[0]; // uu
    tmp[1] = w * std::conj(c[0]) * c[1]; // ud
    tmp[2] = w * std::conj(c[1]) * c[0]; // du
    tmp[3] = w * std::conj(c[1]) * c[1]; // dd
}

// physical magnetization of a normalized spinor: m_i = <c| sigma_i |c>
void physical_m(const cd c[2], double m[3])
{
    m[0] = 2.0 * std::real(std::conj(c[0]) * c[1]);
    m[1] = 2.0 * std::imag(std::conj(c[0]) * c[1]);
    m[2] = std::norm(c[0]) - std::norm(c[1]);
}
} // namespace

TEST(SocMagnetizationRoundtrip, ExtractRecoversPhysicalMagnetization)
{
    // several magnetization directions, all with a nonzero transverse (y) part
    const double dirs[5][3] = {
        {0.0, 1.0, 0.0},                                  // pure +y (the critical case)
        {0.0, -1.0, 0.0},                                 // pure -y (like Mn3Sn atom-1)
        {0.6, 0.8, 0.0},                                  // in-plane 120-deg-like
        {0.36, 0.48, -0.8},                               // general 3D
        {-0.5, 0.5, 0.70710678},                          // general 3D
    };

    // step_trace for a single 2x2 spinor block written contiguously as a 2x2 (col_size=2)
    const int col_size = 2;
    const int step_trace[4] = {0, 1, col_size, col_size + 1};

    for (const auto& mhat : dirs)
    {
        cd c[2];
        spinor_from_direction(mhat, c);

        double m_ref[3];
        physical_m(c, m_ref); // the TRUE magnetization encoded in the state

        cd tmp[4];
        build_DM_block_as_cal_dm_psi(c, 1.0, tmp);

        // 2x2 output buffer (row-major), func writes rho0/x/y/z into step_trace slots at icol=0
        double out[4] = {0, 0, 0, 0};
        elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out);

        const double mx = out[step_trace[1]];
        const double my = out[step_trace[2]];
        const double mz = out[step_trace[3]];

        EXPECT_NEAR(mx, m_ref[0], 1e-10) << "m_x wrong for dir (" << mhat[0] << "," << mhat[1] << "," << mhat[2] << ")";
        EXPECT_NEAR(my, m_ref[1], 1e-10) << "m_y SIGN/VALUE wrong (transverse channel, #7664 regression) for dir ("
                                         << mhat[0] << "," << mhat[1] << "," << mhat[2] << ")";
        EXPECT_NEAR(mz, m_ref[2], 1e-10) << "m_z wrong for dir (" << mhat[0] << "," << mhat[1] << "," << mhat[2] << ")";
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
