#include "mpi.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#define private public
#define protected public
#include "../symmetry_rotation.h"
#undef private
#undef protected
#define DOUBLETHRESHOLD 1e-8

/*

Focused tests for restore_HR_nspin4: four spinor channels of the real-space
EXX H(R), short/long Coulomb channels each restored once, no channel
cross-talk for U = I, correct SU(2) mixing for non-trivial U, and the
antiunitary sigma_y (.)^* sigma_y channel remap (real AND complex inputs).

COVERAGE BOUNDARY (honest): this suite tests the restore_HR_nspin4 helper
directly. The production caller chain Exx_LRI::cal_exx_elec_soc (per-Coulomb-
channel single restore, short/long accumulation) is NOT unit-tested here; it
depends on the full SCF/EXX environment and is covered only indirectly by the
SCF-level cross-feature runs. That call-chain coverage gap is recorded in
MERGE_AUDIT.md.

Conventions (from symmetry_rotation_R.hpp restore_HR_nspin4):

  step 1: each of the 4 irreducible channels is orbitally rotated with
          rotate_atompair_serial(mode='H') -> T1^\dagger H T2 (T real here).
  step 2: SU(2) mixing with U (Su2 = array<complex,4>):
          Hout[a*2+b] = sum_{c,d} conj(U[c*2+a]) * U[d*2+b] * G[c*2+d]
  step 3: antiunitary (isym >= nsym_):
          out[0]=conj(in[3]); out[1]=-conj(in[2]); out[2]=-conj(in[1]);
          out[3]=conj(in[0])

Test matrices (nw = 3; T = rot-90 about z, real Y_l^1 basis):

    T = [ 0 -1  0 ]      A0 = [ 1 2 3 ]   (channel 00)
        [ 1  0  0 ]           [ 4 5 6 ]
        [ 0  0  1 ]           [ 7 8 9 ]

    T^T * A0 * T = [ 5 -4  6 ]
                   [ -2 1 -3 ]
                   [ 8 -7  9 ]

    A1 (channel 01) = 2*A0 ; A2 (channel 10) = 3*A0 ; A3 (channel 11) = 4*A0
    (real; conjugation then acts as identity, keeping the remap check exact)

    Non-trivial SU(2): U = {a,b,-b,a} with a=b=1/sqrt(2) (real, 90 deg):
      Hout[00] = 0.5*(G00 + G01 - G10 - G11)
      Hout[01] = 0.5*(G01 + G00 - G11 - G10)   (same by symmetry)
*/

// mocks
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
RI::Tensor<std::complex<double>> make_T_rot90()
{
    RI::Tensor<std::complex<double>> T({3, 3});
    const double m[3][3] = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T(i, j) = std::complex<double>(m[i][j], 0.0);
    return T;
}

RI::Tensor<std::complex<double>> make_chan(const int scale)
{
    RI::Tensor<std::complex<double>> A({3, 3});
    const double v[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A(i, j) = std::complex<double>(scale * v[i][j], 0.0);
    return A;
}

// complex channel input: real part = scale*v, imaginary part = (scale+1)*0.25
// for the diagonal, (scale+1)*0.5*i off-diagonal - all four channels differ
// in both real and imaginary parts.
RI::Tensor<std::complex<double>> make_chan_complex(const int scale)
{
    RI::Tensor<std::complex<double>> A({3, 3});
    const double v[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A(i, j) = std::complex<double>(scale * v[i][j],
                                           (i == j ? 0.25 : 0.5) * (scale + 1));
    return A;
}

// manual T^T * A * T for the scaled channel
RI::Tensor<std::complex<double>> expected_rotated(const int scale)
{
    RI::Tensor<std::complex<double>> E({3, 3});
    const double v[3][3] = {{5, -4, 6}, {-2, 1, -3}, {8, -7, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            E(i, j) = std::complex<double>(scale * v[i][j], 0.0);
    return E;
}

template <typename Tdata>
bool tensor_close(const RI::Tensor<Tdata>& a, const RI::Tensor<Tdata>& b, double tol)
{
    if (a.shape.size() != b.shape.size())
        return false;
    for (size_t d = 0; d < a.shape.size(); ++d)
        if (a.shape[d] != b.shape[d])
            return false;
    for (size_t i = 0; i < a.shape[0]; ++i)
        for (size_t j = 0; j < a.shape[1]; ++j)
            if (std::abs(a(i, j) - b(i, j)) > tol)
                return false;
    return true;
}
} // namespace

class Nspin4RestoreTest : public testing::Test
{
protected:
    void SetUp() override
    {
        symrot.reduce_Cs_ = true;
        symrot.nsym_ = 1;
        symrot.nanti_ = 1;
        symrot.rotmat_Slm_.resize(2);
        for (int isym = 0; isym < 2; ++isym)
        {
            symrot.rotmat_Slm_[isym].resize(2);
            symrot.rotmat_Slm_[isym][0] = RI::Tensor<std::complex<double>>({1, 1});
            symrot.rotmat_Slm_[isym][0](0, 0) = 1.0;
            symrot.rotmat_Slm_[isym][1] = make_T_rot90();
        }
        // identity SU(2) on both operations by default
        symrot.spin_U_.resize(2);
        for (int isym = 0; isym < 2; ++isym)
            symrot.spin_U_[isym] = ModuleSymmetry::SpinRotation::Su2{1.0, 0.0, 0.0, 1.0};

        atoms[0].nw = 3;
        atoms[1].nw = 3;
        atoms[0].label = "H";
        atoms[1].label = "H";
        // one l=1 shell: iw2l = {1,1,1} so set_rotation_matrix builds the 3x3 block
        atoms[0].iw2l = {1, 1, 1};
        atoms[1].iw2l = {1, 1, 1};

        st.iat2it = new int[2]{0, 0};
    }
    // NOTE: Statistics' destructor deletes iat2it itself; do NOT delete here.

    std::array<std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<std::complex<double>>>>, 4>
    make_irr_input()
    {
        std::array<std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<std::complex<double>>>>, 4> in;
        for (int is = 0; is < 4; ++is)
            in[is][0][{0, {0, 0, 0}}] = make_chan(is + 1);
        return in;
    }

    ModuleSymmetry::Symmetry_rotation symrot;
    Atom atoms[2];
    Statistics st;
};

TEST_F(Nspin4RestoreTest, FourChannelsIndependentForIdentitySU2)
{
    // star with a single unitary member: 4 channels restored independently,
    // each equal to T^T * A_channel * T - no spin mixing, no cross-talk.
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {{0, {{0, 0}, {0, 0, 0}}}};

    auto in = make_irr_input();
    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);

    for (int is = 0; is < 4; ++is)
    {
        EXPECT_EQ(out[is].size(), 1u);
        EXPECT_TRUE(tensor_close(out[is][0][{0, {0, 0, 0}}], expected_rotated(is + 1), DOUBLETHRESHOLD));
    }
}

TEST_F(Nspin4RestoreTest, OffDiagonalChannelsNonzero)
{
    // all four channels (including 01 and 10) are nonzero and distinct in
    // both the input and the restored output.
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {{0, {{0, 0}, {0, 0, 0}}}};
    auto in = make_irr_input();
    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);
    for (int is = 0; is < 4; ++is)
    {
        const auto& H = out[is][0][{0, {0, 0, 0}}];
        bool nonzero = false;
        for (size_t i = 0; i < 3 && !nonzero; ++i)
            for (size_t j = 0; j < 3 && !nonzero; ++j)
                nonzero = std::abs(H(i, j)) > 1e-6;
        EXPECT_TRUE(nonzero) << "channel " << is << " became zero";
    }
}

TEST_F(Nspin4RestoreTest, NonTrivialSU2Mixing)
{
    // U = {a,b,-b,a} with a=b=1/sqrt(2). From
    // Hout[a*2+b] = sum_{c,d} conj(U[c*2+a]) * U[d*2+b] * G[c*2+d]:
    //   Hout[00] = 0.5*(G00 - G01 - G10 + G11)
    //   Hout[01] = 0.5*(G00 + G01 - G10 - G11)
    symrot.spin_U_[0] = ModuleSymmetry::SpinRotation::Su2{M_SQRT1_2, M_SQRT1_2, -M_SQRT1_2, M_SQRT1_2};
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {{0, {{0, 0}, {0, 0, 0}}}};

    auto in = make_irr_input();
    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);

    // G[is] = T^T * (scale is+1) * A0 * T = (is+1) * expected_rotated(1)
    RI::Tensor<std::complex<double>> G0 = expected_rotated(1);
    RI::Tensor<std::complex<double>> G1 = expected_rotated(2);
    RI::Tensor<std::complex<double>> G2 = expected_rotated(3);
    RI::Tensor<std::complex<double>> G3 = expected_rotated(4);
    RI::Tensor<std::complex<double>> E00({3, 3});
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            E00(i, j) = 0.5 * (G0(i, j) - G1(i, j) - G2(i, j) + G3(i, j));
    EXPECT_TRUE(tensor_close(out[0][0][{0, {0, 0, 0}}], E00, DOUBLETHRESHOLD));
    RI::Tensor<std::complex<double>> E01({3, 3});
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            E01(i, j) = 0.5 * (G0(i, j) + G1(i, j) - G2(i, j) - G3(i, j));
    EXPECT_TRUE(tensor_close(out[1][0][{0, {0, 0, 0}}], E01, DOUBLETHRESHOLD));
}

TEST_F(Nspin4RestoreTest, AntiunitarySigmaYChannelRemap)
{
    // star member isym=1 is antiunitary (>= nsym_): after the (identity)
    // SU(2) step, channels remap as out[0]=conj(in[3]), out[1]=-conj(in[2]),
    // out[2]=-conj(in[1]), out[3]=conj(in[0]).
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {
        {0, {{0, 0}, {0, 0, 0}}},
        {1, {{1, 1}, {1, 0, 0}}},
    };

    auto in = make_irr_input();
    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);

    // member 0 (unitary): independent channels
    for (int is = 0; is < 4; ++is)
        EXPECT_TRUE(tensor_close(out[is][0][{0, {0, 0, 0}}], expected_rotated(is + 1), DOUBLETHRESHOLD));
    // member 1 (antiunitary): sigma_y K remap, real input -> conjugation identity
    RI::Tensor<std::complex<double>> neg3 = expected_rotated(3);
    RI::Tensor<std::complex<double>> neg2 = expected_rotated(2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            neg3(i, j) = -neg3(i, j);
            neg2(i, j) = -neg2(i, j);
        }
    EXPECT_TRUE(tensor_close(out[0][1][{1, {1, 0, 0}}], expected_rotated(4), DOUBLETHRESHOLD)); // 00 <- +conj(11)
    EXPECT_TRUE(tensor_close(out[1][1][{1, {1, 0, 0}}], neg3, DOUBLETHRESHOLD));               // 01 <- -conj(10)
    EXPECT_TRUE(tensor_close(out[2][1][{1, {1, 0, 0}}], neg2, DOUBLETHRESHOLD));               // 10 <- -conj(01)
    EXPECT_TRUE(tensor_close(out[3][1][{1, {1, 0, 0}}], expected_rotated(1), DOUBLETHRESHOLD)); // 11 <- +conj(00)
}

TEST_F(Nspin4RestoreTest, AntiunitarySigmaYComplexInputs)
{
    // complex inputs: the antiunitary branch must conjugate elementwise AND
    // apply the channel remap with sign flips. With U = I:
    //   out[0] = conj(rot(A3)), out[1] = -conj(rot(A2)),
    //   out[2] = -conj(rot(A1)), out[3] = conj(rot(A0))
    // where rot(Ak) = T^T * A_k * T is the orbitally rotated channel.
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {
        {0, {{0, 0}, {0, 0, 0}}},
        {1, {{1, 1}, {1, 0, 0}}},
    };

    std::array<std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<std::complex<double>>>>, 4> in;
    for (int is = 0; is < 4; ++is)
        in[is][0][{0, {0, 0, 0}}] = make_chan_complex(is + 1);

    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);

    // unitary member (isym=0): channels pass through unchanged (no conj)
    for (int is = 0; is < 4; ++is)
    {
        // rot(Ak) computed manually: T^T A T with T = rot-90; for A_k with
        // entries a_ij: (T^T A T)_00 = a_11, (T^T A T)_01 = -a_10,
        // (T^T A T)_02 = a_12, (T^T A T)_10 = -a_01, (T^T A T)_11 = a_00,
        // (T^T A T)_12 = -a_02, (T^T A T)_20 = a_21, (T^T A T)_21 = -a_20,
        // (T^T A T)_22 = a_22
        RI::Tensor<std::complex<double>> E({3, 3});
        const auto& A = in[is][0][{0, {0, 0, 0}}];
        E(0, 0) = A(1, 1); E(0, 1) = -A(1, 0); E(0, 2) = A(1, 2);
        E(1, 0) = -A(0, 1); E(1, 1) = A(0, 0); E(1, 2) = -A(0, 2);
        E(2, 0) = A(2, 1); E(2, 1) = -A(2, 0); E(2, 2) = A(2, 2);
        const bool ok_u = tensor_close(out[is][0][{0, {0, 0, 0}}], E, DOUBLETHRESHOLD);
        EXPECT_TRUE(ok_u) << "unitary member channel " << is;
    }
    // antiunitary member (isym=1): conjugate + remap + sign
    for (int is = 0; is < 4; ++is)
    {
        // out[0] <- +conj(rot(A3)); out[1] <- -conj(rot(A2));
        // out[2] <- -conj(rot(A1)); out[3] <- +conj(rot(A0))
        const int src[4] = {3, 2, 1, 0};
        const bool neg[4] = {false, true, true, false};
        // rotate channel src[is] (input index src[is] -> make_chan_complex(src[is]+1))
        const auto& Asrc = in[src[is]][0][{0, {0, 0, 0}}];
        RI::Tensor<std::complex<double>> rotc({3, 3});
        rotc(0, 0) = Asrc(1, 1); rotc(0, 1) = -Asrc(1, 0); rotc(0, 2) = Asrc(1, 2);
        rotc(1, 0) = -Asrc(0, 1); rotc(1, 1) = Asrc(0, 0); rotc(1, 2) = -Asrc(0, 2);
        rotc(2, 0) = Asrc(2, 1); rotc(2, 1) = -Asrc(2, 0); rotc(2, 2) = Asrc(2, 2);
        RI::Tensor<std::complex<double>> E({3, 3});
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                E(i, j) = neg[is] ? -std::conj(rotc(i, j)) : std::conj(rotc(i, j));
        const bool ok_a = tensor_close(out[is][1][{1, {1, 0, 0}}], E, DOUBLETHRESHOLD);
        EXPECT_TRUE(ok_a) << "antiunitary member channel " << is;
    }
    // explicit sanity: out[3] <- +conj(rot(A0)) (src[3]=0, neg[3]=false).
    // rot(A0)(0,0) = A0(1,1) = 5 + 0.5i  ->  conj -> 5 - 0.5i
    const auto& out3 = out[3][1][{1, {1, 0, 0}}];
    EXPECT_NEAR(out3(0, 0).real(), 5.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(out3(0, 0).imag(), -0.5, DOUBLETHRESHOLD);
}

TEST_F(Nspin4RestoreTest, ShortAndLongRestoreOnceEach)
{
    // short and long Coulomb channels are restored by two separate calls;
    // each call restores all four spin channels exactly once. Verify by
    // running two independent calls and checking per-call output size.
    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {{0, {{0, 0}, {0, 0, 0}}}};
    auto in_short = make_irr_input();
    auto in_long = make_irr_input();
    // long channel carries different values to prove no cross-talk between calls
    for (int is = 0; is < 4; ++is)
        in_long[is][0][{0, {0, 0, 0}}] = make_chan(10 * (is + 1));

    ModuleSymmetry::Symmetry symm;
    auto out_short = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in_short);
    auto out_long = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in_long);

    for (int is = 0; is < 4; ++is)
    {
        EXPECT_TRUE(tensor_close(out_short[is][0][{0, {0, 0, 0}}], expected_rotated(is + 1), DOUBLETHRESHOLD));
        EXPECT_TRUE(tensor_close(out_long[is][0][{0, {0, 0, 0}}], expected_rotated(10 * (is + 1)), DOUBLETHRESHOLD));
    }
}

TEST_F(Nspin4RestoreTest, HermiticityPreservedForUnitaryMember)
{
    // Hermitian input channel 00 -> Hermitian output under the unitary member
    // (T real orthogonal, U=I: conjugation-free path).
    RI::Tensor<std::complex<double>> H({3, 3});
    H(0, 0) = std::complex<double>(2, 0); H(0, 1) = std::complex<double>(1, 1); H(0, 2) = std::complex<double>(3, -2);
    H(1, 0) = std::conj(H(0, 1)); H(1, 1) = std::complex<double>(4, 0); H(1, 2) = std::complex<double>(0.5, 0.25);
    H(2, 0) = std::conj(H(0, 2)); H(2, 1) = std::conj(H(1, 2)); H(2, 2) = std::complex<double>(6, 0);

    symrot.irs_.sector_stars_[{{0, 0}, {0, 0, 0}}] = {{0, {{0, 0}, {0, 0, 0}}}};
    std::array<std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<std::complex<double>>>>, 4> in;
    in[0][0][{0, {0, 0, 0}}] = H;

    ModuleSymmetry::Symmetry symm;
    auto out = symrot.restore_HR_nspin4(symm, atoms, st, 'H', in);

    const auto& Hout = out[0][0][{0, {0, 0, 0}}];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(std::abs(Hout(i, j) - std::conj(Hout(j, i))), 0.0, DOUBLETHRESHOLD);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
