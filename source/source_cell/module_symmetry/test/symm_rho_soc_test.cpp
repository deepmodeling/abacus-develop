#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <array>

#include "../symmetry.h"
#include "../symm_rot_spin.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/density_matrix.h" // real func_xyz_to_updown

/************************************************
 *  unit test of Symmetry::rhog_symmetry_nspin4
 *  (nspin=4 / SOC reciprocal-space spin-density symmetrization)
 *
 *  The operator is driven with a MANUALLY built symmetry group (no analy_sys),
 *  so the group, grid and spin rotations are fully controlled. We use the proper
 *  point group D_4 = {E, C4z, C2z, C4z^3, C2x, C2y, C2[110], C2[1-10]} on a cubic
 *  lattice (a=1 => gmatc = kgmatrix = gmatrix = R^T). D_4 is NON-ABELIAN, so the
 *  test is sensitive to spin-rotation representation/handedness bugs that an
 *  abelian group (e.g. C_6h) would hide.
 *
 *  Checks:
 *   - Idempotence: symmetrizing an already-symmetric density is a no-op.
 *   - Invariance:  the symmetrized density satisfies m(R_g G) = W(g) m(G) for
 *                  every group operation g, verified by an independent oracle.
***********************************************/

// mock the unused constructors pulled in by linking the symmetry library
pseudo::pseudo() {}
pseudo::~pseudo() {}
Atom::Atom() {}
Atom::~Atom() {}
Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

namespace
{
constexpr int N = 4;                // grid dimension (even, so -i mod N stays on-grid)
constexpr int NXYZ = N * N * N;
constexpr double TOL = 1e-10;

// the 8 proper rotations of D_4, as textbook column-vector matrices R (r'=R r)
const std::array<std::array<std::array<int, 3>, 3>, 8> Rcol = {{
    {{{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}}},   // E
    {{{ 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}}},   // C4z
    {{{-1, 0, 0}, { 0,-1, 0}, { 0, 0, 1}}},   // C2z
    {{{ 0, 1, 0}, {-1, 0, 0}, { 0, 0, 1}}},   // C4z^3
    {{{ 1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}}},   // C2x
    {{{-1, 0, 0}, { 0, 1, 0}, { 0, 0,-1}}},   // C2y
    {{{ 0, 1, 0}, { 1, 0, 0}, { 0, 0,-1}}},   // C2[110]
    {{{ 0,-1, 0}, {-1, 0, 0}, { 0, 0,-1}}},   // C2[1-10]
}};

// ABACUS stores the cartesian rotation in the row-vector convention: gmatc = R^T.
ModuleBase::Matrix3 gmatc_of(int g)
{
    const auto& R = Rcol[g];
    return ModuleBase::Matrix3(R[0][0], R[1][0], R[2][0],
                               R[0][1], R[1][1], R[2][1],
                               R[0][2], R[1][2], R[2][2]);
}

// mirror of the internal rotate_recip: G' index components from kgmatrix (=gmatc)
void rotate_index(const ModuleBase::Matrix3& g, int i, int j, int k, int& ii, int& jj, int& kk)
{
    ii = int(g.e11 * i + g.e21 * j + g.e31 * k); if (ii < 0) ii += 10 * N; ii %= N;
    jj = int(g.e12 * i + g.e22 * j + g.e32 * k); if (jj < 0) jj += 10 * N; jj %= N;
    kk = int(g.e13 * i + g.e23 * j + g.e33 * k); if (kk < 0) kk += 10 * N; kk %= N;
}

// build a Symmetry object carrying the D_4 group on a cubic grid
void build_group(ModuleSymmetry::Symmetry& symm, std::vector<ModuleBase::Matrix3>& wspin)
{
    symm.epsilon = 1e-6;
    symm.nrot = 8;
    symm.nrotk = 8;
    symm.ncell = 1;
    symm.ptrans = {ModuleBase::Vector3<double>(0.0, 0.0, 0.0)};
    ModuleSymmetry::Symmetry::pricell_loop = false;
    wspin.resize(8);
    for (int g = 0; g < 8; ++g)
    {
        const ModuleBase::Matrix3 gc = gmatc_of(g);
        symm.gmatrix[g] = gc;      // cubic a=1: direct == cartesian
        symm.kgmatrix[g] = gc;     // orthogonal rotation: reciprocal == direct
        symm.gtrans[g] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        wspin[g] = ModuleSymmetry::SpinRotation::spin_so3(gc);
    }
}

// a fixed, non-symmetric complex spin density on the grid
void fill_density(std::vector<std::complex<double>>& x,
                  std::vector<std::complex<double>>& y,
                  std::vector<std::complex<double>>& z)
{
    for (int idx = 0; idx < NXYZ; ++idx)
    {
        x[idx] = std::complex<double>(0.3 * idx - 1.0, 0.7 * ((idx * 13) % 5) - 1.5);
        y[idx] = std::complex<double>(-0.5 * ((idx * 7) % 4) + 0.9, 0.2 * idx - 2.0);
        z[idx] = std::complex<double>(0.11 * ((idx * 3) % 6), -0.4 * ((idx * 5) % 7) + 1.0);
    }
}
} // namespace

TEST(RhogSymmetrySoc, Idempotence)
{
    ModuleSymmetry::Symmetry symm;
    std::vector<ModuleBase::Matrix3> wspin;
    build_group(symm, wspin);

    std::vector<int> ixyz2ipw(NXYZ);
    for (int i = 0; i < NXYZ; ++i) { ixyz2ipw[i] = i; } // every FFT point is a plane wave

    std::vector<std::complex<double>> x(NXYZ), y(NXYZ), z(NXYZ);
    fill_density(x, y, z);

    symm.rhog_symmetry_nspin4(x.data(), y.data(), z.data(), wspin.data(), ixyz2ipw.data(), N, N, N, N, N, N, nullptr, nullptr, nullptr, -1);
    std::vector<std::complex<double>> x1 = x, y1 = y, z1 = z;
    symm.rhog_symmetry_nspin4(x.data(), y.data(), z.data(), wspin.data(), ixyz2ipw.data(), N, N, N, N, N, N, nullptr, nullptr, nullptr, -1);

    for (int i = 0; i < NXYZ; ++i)
    {
        EXPECT_NEAR(x[i].real(), x1[i].real(), TOL); EXPECT_NEAR(x[i].imag(), x1[i].imag(), TOL);
        EXPECT_NEAR(y[i].real(), y1[i].real(), TOL); EXPECT_NEAR(y[i].imag(), y1[i].imag(), TOL);
        EXPECT_NEAR(z[i].real(), z1[i].real(), TOL); EXPECT_NEAR(z[i].imag(), z1[i].imag(), TOL);
    }
}

TEST(RhogSymmetrySoc, GroupInvariance)
{
    ModuleSymmetry::Symmetry symm;
    std::vector<ModuleBase::Matrix3> wspin;
    build_group(symm, wspin);

    std::vector<int> ixyz2ipw(NXYZ);
    for (int i = 0; i < NXYZ; ++i) { ixyz2ipw[i] = i; }

    std::vector<std::complex<double>> x(NXYZ), y(NXYZ), z(NXYZ);
    fill_density(x, y, z);
    symm.rhog_symmetry_nspin4(x.data(), y.data(), z.data(), wspin.data(), ixyz2ipw.data(), N, N, N, N, N, N, nullptr, nullptr, nullptr, -1);

    // non-triviality guard: the symmetrized density must not be all-zero, otherwise
    // invariance would hold trivially and the test would be meaningless.
    double maxabs = 0.0;
    for (int i = 0; i < NXYZ; ++i)
    {
        maxabs = std::max(maxabs, std::abs(x[i]));
        maxabs = std::max(maxabs, std::abs(y[i]));
        maxabs = std::max(maxabs, std::abs(z[i]));
    }
    EXPECT_GT(maxabs, 0.1);

    // independent oracle: the symmetrized density must obey m(R_g G) = W(g) m(G) for all g, G.
    for (int g = 0; g < 8; ++g)
    {
        const ModuleBase::Matrix3& W = wspin[g];
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    const int idx = (i * N + j) * N + k;
                    int ii, jj, kk;
                    rotate_index(symm.kgmatrix[g], i, j, k, ii, jj, kk);
                    const int idx2 = (ii * N + jj) * N + kk;
                    const std::complex<double> ex = W.e11 * x[idx] + W.e12 * y[idx] + W.e13 * z[idx];
                    const std::complex<double> ey = W.e21 * x[idx] + W.e22 * y[idx] + W.e23 * z[idx];
                    const std::complex<double> ez = W.e31 * x[idx] + W.e32 * y[idx] + W.e33 * z[idx];
                    EXPECT_NEAR(x[idx2].real(), ex.real(), TOL) << "g=" << g << " idx=" << idx;
                    EXPECT_NEAR(x[idx2].imag(), ex.imag(), TOL) << "g=" << g << " idx=" << idx;
                    EXPECT_NEAR(y[idx2].real(), ey.real(), TOL) << "g=" << g << " idx=" << idx;
                    EXPECT_NEAR(y[idx2].imag(), ey.imag(), TOL) << "g=" << g << " idx=" << idx;
                    EXPECT_NEAR(z[idx2].real(), ez.real(), TOL) << "g=" << g << " idx=" << idx;
                    EXPECT_NEAR(z[idx2].imag(), ez.imag(), TOL) << "g=" << g << " idx=" << idx;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Coupling test (nonzero m_y): the spin-density rotation W=spin_so3 used by psymmg_soc for the
// grid symmetrization MUST agree with the SU(2) rotation of the physical spinor state followed by
// the REAL func_xyz_to_updown extraction (which reads the conj-first stored DM, DM=conj(P), and
// uses the bare +Im(ud)-Im(du)). This test now calls the actual func_xyz_to_updown rather than a
// local re-implementation, so the grid-rotation and DM-extraction conventions cannot drift apart
// silently (it fails on the #7664 m_y flip). The self-referential GroupInvariance test above
// cannot catch this because it uses the same wspin as its own oracle.
// ---------------------------------------------------------------------------
namespace
{
using cd = std::complex<double>;
// PHYSICAL spinor block P = r0*I + m.sigma  (sigma_y = [[0,-i],[i,0]]); layout {uu,ud,du,dd}
ModuleSymmetry::SpinRotation::Su2 block_from_pauli(double r0, double mx, double my, double mz)
{
    return {cd(r0 + mz, 0.0), cd(mx, -my), cd(mx, my), cd(r0 - mz, 0.0)};
}
// The runtime stores the DM conj-first (DM = conj(P), cal_dm_psi); this is what func_xyz_to_updown
// actually consumes. Given a physical block P, the stored block is its element-wise conjugate.
ModuleSymmetry::SpinRotation::Su2 stored_dm_from_phys(const ModuleSymmetry::SpinRotation::Su2& P)
{
    return {std::conj(P[0]), std::conj(P[1]), std::conj(P[2]), std::conj(P[3])};
}
// call the REAL func_xyz_to_updown on a 2x2 stored-DM block; return (m_x, m_y, m_z)
ModuleBase::Vector3<double> real_extract(const ModuleSymmetry::SpinRotation::Su2& Dstored)
{
    const cd tmp[4] = {Dstored[0], Dstored[1], Dstored[2], Dstored[3]}; // {uu,ud,du,dd}
    const int col_size = 2;
    const int step_trace[4] = {0, 1, col_size, col_size + 1};
    double out[4] = {0.0, 0.0, 0.0, 0.0}; // rho0/x/y/z written at icol=0
    elecstate::DensityMatrix_Tools::func_xyz_to_updown<double>(tmp, 0, step_trace, out);
    return ModuleBase::Vector3<double>(out[step_trace[1]], out[step_trace[2]], out[step_trace[3]]);
}
} // namespace

TEST(RhogSymmetrySoc, SpinConventionCoupling)
{
    // representative magnetizations, all with a nonzero y-component
    const double mtest[4][3] = {{0.4, 0.7, -0.5}, {0.0, 1.0, 0.0}, {-0.3, 0.6, 0.9}, {1.0, -0.8, 0.2}};

    for (int g = 0; g < 8; ++g)
    {
        const ModuleBase::Matrix3 gc = gmatc_of(g);
        const ModuleSymmetry::SpinRotation::Su2 U = ModuleSymmetry::SpinRotation::so3_to_su2(gc);
        const ModuleBase::Matrix3 Wgrid = ModuleSymmetry::SpinRotation::spin_so3(gc);
        const ModuleBase::Matrix3 Wpauli = ModuleSymmetry::SpinRotation::pauli_rotation_matrix(U);

        // (1) the geometric grid rotation and the SU(2)-induced Pauli rotation must coincide
        EXPECT_NEAR(Wgrid.e11, Wpauli.e11, TOL) << "g=" << g; EXPECT_NEAR(Wgrid.e12, Wpauli.e12, TOL) << "g=" << g;
        EXPECT_NEAR(Wgrid.e13, Wpauli.e13, TOL) << "g=" << g; EXPECT_NEAR(Wgrid.e21, Wpauli.e21, TOL) << "g=" << g;
        EXPECT_NEAR(Wgrid.e22, Wpauli.e22, TOL) << "g=" << g; EXPECT_NEAR(Wgrid.e23, Wpauli.e23, TOL) << "g=" << g;
        EXPECT_NEAR(Wgrid.e31, Wpauli.e31, TOL) << "g=" << g; EXPECT_NEAR(Wgrid.e32, Wpauli.e32, TOL) << "g=" << g;
        EXPECT_NEAR(Wgrid.e33, Wpauli.e33, TOL) << "g=" << g;

        // (2) End-to-end with the REAL func_xyz_to_updown, exactly the runtime data flow:
        //   physical block P(m) --conj--> stored DM (conj-first) --func_xyz_to_updown--> grid m.
        //   Rotate the PHYSICAL block by the spinor SU(2) U (U P U^dagger, i.e. the physical state
        //   rotation), conj to the stored block, extract again -> m'. psymmg_soc rotates the grid
        //   components with Wgrid=spin_so3, so we must have  m' == Wgrid * m.  This catches any
        //   mismatch (e.g. the #7664 m_y flip) between func_xyz_to_updown and spin_so3.
        for (const auto& m : mtest)
        {
            const ModuleSymmetry::SpinRotation::Su2 P  = block_from_pauli(2.0, m[0], m[1], m[2]);
            const ModuleSymmetry::SpinRotation::Su2 Pp = ModuleSymmetry::SpinRotation::rotate_spin_block(P, U);

            const ModuleBase::Vector3<double> mF  = real_extract(stored_dm_from_phys(P));
            const ModuleBase::Vector3<double> mFp = real_extract(stored_dm_from_phys(Pp));

            // (2a) extraction recovers the physical magnetization (block carries 2*m)
            EXPECT_NEAR(mF.x, 2.0 * m[0], TOL) << "g=" << g;
            EXPECT_NEAR(mF.y, 2.0 * m[1], TOL) << "g=" << g << " (m_y extraction)";
            EXPECT_NEAR(mF.z, 2.0 * m[2], TOL) << "g=" << g;

            // (2b) grid rotation spin_so3 agrees with the SU(2) block rotation + real extraction
            const ModuleBase::Vector3<double> mrot = Wgrid * mF;
            EXPECT_NEAR(mFp.x, mrot.x, TOL) << "g=" << g;
            EXPECT_NEAR(mFp.y, mrot.y, TOL) << "g=" << g << " (y-channel handedness)";
            EXPECT_NEAR(mFp.z, mrot.z, TOL) << "g=" << g;
        }
    }
}

int main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
