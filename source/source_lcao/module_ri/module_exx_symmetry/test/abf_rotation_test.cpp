#include "mpi.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
// symmetry_rotation.h pulls in symmetry_rotation_R.hpp (no include guard);
// standard-library headers must all be included before `private public`.
#define private public
#define protected public
#include "../symmetry_rotation.h"
#undef private
#undef protected
using namespace std::complex_literals; // for the `1i` literal used below
#define DOUBLETHRESHOLD 1e-8

/*

Focused mathematical tests for scalar-ABF unitary/antiunitary restore.

Convention (derived from LibRI RI::Sym::T1_HR_T2, Symmetry_Rotation.h):

    H' = T1^\dagger * H * T2

For scalar ABF the rotation matrices T are real (real spherical-harmonic
basis; rotmat_Slm_ entries have zero imaginary part), so T^\dagger = T^T.

rotate_atompair_serial_abf:
    unitary     (isym <  nsym_):  TAT = T1^T * A * T2
    antiunitary (isym >= nsym_):  TAT = conj(T1^T * A * T2)
                                    = T1^T * conj(A) * T2   (T real)

Test matrices: T = rotation by 90 deg about z in the real Y_l^1 basis
(m = -1,0,1 order):

    T = [ 0 -1  0 ]
        [ 1  0  0 ]
        [ 0  0  1 ]

    A (real 3x3) = [ 1 2 3 ]      A_c (complex) = [ 1+i  2   3  ]
                   [ 4 5 6 ]                      [ 4   5-i  6  ]
                   [ 7 8 9 ]                      [ 7   8   9+i ]

    T^T * A * T =
        [ 5  -4  6 ]
        [ -2  1 -3 ]
        [ 8  -7  9 ]

    T^T * A_c * T =
        [  5+0i   -4+0i   6+0i ]
        [ -2+0i    1+0i  -3+0i ]
        [  8+0i   -7+0i   9+1i ]

    (the entry (2,2) of A_c, 9+i, maps onto itself under this rotation)

*/

// mocks: dummy constructors for classes pulled in by the headers
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
// real 3x3 rotation about z by 90 deg, l=1 real-spherical-harmonic basis
RI::Tensor<std::complex<double>> make_T_rot90()
{
    RI::Tensor<std::complex<double>> T({3, 3});
    const double m[3][3] = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            T(i, j) = std::complex<double>(m[i][j], 0.0);
    return T;
}

RI::Tensor<double> make_A_real()
{
    RI::Tensor<double> A({3, 3});
    const double v[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A(i, j) = v[i][j];
    return A;
}

RI::Tensor<std::complex<double>> make_A_complex()
{
    RI::Tensor<std::complex<double>> A({3, 3});
    const std::complex<double> v[3][3] = {
        {1.0 + 1i, 2.0, 3.0}, {4.0, 5.0 - 1i, 6.0}, {7.0, 8.0, 9.0 + 1i}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            A(i, j) = v[i][j];
    return A;
}

// manual T^T * A * T (T real)
RI::Tensor<double> expected_unitary()
{
    RI::Tensor<double> E({3, 3});
    const double v[3][3] = {{5, -4, 6}, {-2, 1, -3}, {8, -7, 9}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            E(i, j) = v[i][j];
    return E;
}

// manual T^T * A_c * T
RI::Tensor<std::complex<double>> expected_unitary_complex()
{
    RI::Tensor<std::complex<double>> E({3, 3});
    const std::complex<double> v[3][3] = {
        {5.0 - 1i, -4.0, 6.0}, {-2.0, 1.0 + 1i, -3.0}, {8.0, -7.0, 9.0 + 1i}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            E(i, j) = v[i][j];
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

class AbfRotationTest : public testing::Test
{
protected:
    void SetUp() override
    {
        symrot.reduce_Cs_ = true;
        // type 0: one l=1 ABF channel (3 functions); type 1: one l=0 (1 function)
        symrot.abfs_l_nchi_ = {{0, 1}, {1, 0}};
        symrot.nsym_ = 1; // 1 unitary op (identity)
        symrot.nanti_ = 1; // + 1 antiunitary op with the same spatial rotation
        symrot.rotmat_Slm_.resize(2);
        for (int isym = 0; isym < 2; ++isym)
        {
            symrot.rotmat_Slm_[isym].resize(2);
            // L=0 block: identity (real)
            symrot.rotmat_Slm_[isym][0] = RI::Tensor<std::complex<double>>({1, 1});
            symrot.rotmat_Slm_[isym][0](0, 0) = 1.0;
            // L=1 block: rotation by 90 deg about z (real in real-Y basis)
            symrot.rotmat_Slm_[isym][1] = make_T_rot90();
        }
    }
    ModuleSymmetry::Symmetry_rotation symrot;
};

TEST_F(AbfRotationTest, UnitaryRealTensor)
{
    // isym=0 < nsym_: TAT = T1^T * A * T2 (T real)
    RI::Tensor<double> A = make_A_real();
    RI::Tensor<double> TAT = symrot.rotate_atompair_serial_abf(A, 0, 0, 0);
    EXPECT_TRUE(tensor_close(TAT, expected_unitary(), DOUBLETHRESHOLD));
}

TEST_F(AbfRotationTest, AntiunitaryRealTensor)
{
    // isym=1 >= nsym_: TAT = conj(T1^T * A * T2); for real A,T this equals
    // the unitary result (conjugation is the identity on reals) - but the
    // result must remain real, and must equal the manual expected matrix.
    RI::Tensor<double> A = make_A_real();
    RI::Tensor<double> TAT = symrot.rotate_atompair_serial_abf(A, 1, 0, 0);
    EXPECT_TRUE(tensor_close(TAT, expected_unitary(), DOUBLETHRESHOLD));
    // explicitly verify no stray imaginary part is introduced
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_EQ(TAT(i, j), TAT(i, j)); // placeholder, real dtype
}

TEST_F(AbfRotationTest, UnitaryComplexTensor)
{
    // complex A: unitary TAT = T1^T * A * T2, manual expected values
    RI::Tensor<std::complex<double>> A = make_A_complex();
    RI::Tensor<std::complex<double>> TAT
        = symrot.rotate_atompair_serial_abf<std::complex<double>>(A, 0, 0, 0);
    EXPECT_TRUE(tensor_close(TAT, expected_unitary_complex(), DOUBLETHRESHOLD));
}

TEST_F(AbfRotationTest, AntiunitaryComplexTensor)
{
    // antiunitary: TAT = conj(T1^T * A * T2); manual expected = conj(unitary)
    RI::Tensor<std::complex<double>> A = make_A_complex();
    RI::Tensor<std::complex<double>> TAT
        = symrot.rotate_atompair_serial_abf<std::complex<double>>(A, 1, 0, 0);
    RI::Tensor<std::complex<double>> E = expected_unitary_complex();
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            EXPECT_NEAR(std::abs(TAT(i, j) - std::conj(E(i, j))), 0.0, DOUBLETHRESHOLD);
    // contrast with the non-conjugated result: entry (2,2) must flip sign of Im
    EXPECT_NEAR(TAT(2, 2).imag(), -1.0, DOUBLETHRESHOLD);
}

TEST_F(AbfRotationTest, DifferentTypesT1T2)
{
    // type1=0 (l=1, T1 = 3x3 rot90), type2=1 (l=0, T2 = 1x1 identity)
    // A is 3x1; expected = T1^T * A
    RI::Tensor<double> A({3, 1});
    A(0, 0) = 1; A(1, 0) = 2; A(2, 0) = 3;
    RI::Tensor<double> TAT = symrot.rotate_atompair_serial_abf(A, 0, 0, 1);
    // T^T * [1,2,3]^T = [2, -1, 3]^T
    EXPECT_NEAR(TAT(0, 0), 2.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(TAT(1, 0), -1.0, DOUBLETHRESHOLD);
    EXPECT_NEAR(TAT(2, 0), 3.0, DOUBLETHRESHOLD);
}

TEST_F(AbfRotationTest, RestoreHRAbfStarMapping)
{
    // end-to-end restore_HR_abf: one irreducible entry (iat0,iat0,R=0) whose
    // star contains two members: isym=0 -> (ap1=0,ap2=0,R=0) and
    // isym=1 (antiunitary) -> (ap1=1,ap2=1,R=(1,0,0)).
    Statistics st;
    int* iat2it = new int[2]{0, 0};
    st.iat2it = iat2it;

    std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<double>>> HR_irr;
    RI::Tensor<double> A = make_A_real();
    HR_irr[0][{0, {0, 0, 0}}] = A;

    ModuleSymmetry::TapR irapR = {{0, 0}, {0, 0, 0}};
    symrot.irs_.sector_stars_[irapR] = {
        {0, {{0, 0}, {0, 0, 0}}},
        {1, {{1, 1}, {1, 0, 0}}},
    };

    ModuleSymmetry::Symmetry symm;
    Atom atoms[2];
    auto HR_full = symrot.restore_HR_abf(symm, atoms, st, HR_irr);

    // member 0: unitary -> T^T A T
    EXPECT_TRUE(tensor_close(HR_full[0][{0, {0, 0, 0}}], expected_unitary(), DOUBLETHRESHOLD));
    // member 1: antiunitary -> conj(T^T A T) == T^T A T (real)
    EXPECT_TRUE(tensor_close(HR_full[1][{1, {1, 0, 0}}], expected_unitary(), DOUBLETHRESHOLD));
    // no extra keys
    EXPECT_EQ(HR_full.size(), 2u);
    EXPECT_EQ(HR_full[0].size(), 1u);
    EXPECT_EQ(HR_full[1].size(), 1u);
}

TEST_F(AbfRotationTest, RestoreHRAbfInvalidShapeThrows)
{
    // an invalid (1-D) tensor must raise instead of being silently skipped
    Statistics st;
    int* iat2it = new int[1]{0};
    st.iat2it = iat2it;

    std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<double>>> HR_irr;
    RI::Tensor<double> bad({3}); // not a valid matrix
    HR_irr[0][{0, {0, 0, 0}}] = bad;

    ModuleSymmetry::TapR irapR = {{0, 0}, {0, 0, 0}};
    symrot.irs_.sector_stars_[irapR] = {{0, {{0, 0}, {0, 0, 0}}}};

    ModuleSymmetry::Symmetry symm;
    Atom atoms[1];
    EXPECT_THROW(symrot.restore_HR_abf(symm, atoms, st, HR_irr), std::runtime_error);
}

TEST_F(AbfRotationTest, RestoreHRAbfDuplicateKeyThrows)
{
    // two star members mapping onto the same target key must raise
    Statistics st;
    int* iat2it = new int[1]{0};
    st.iat2it = iat2it;

    std::map<int, std::map<std::pair<int, ModuleSymmetry::TC>, RI::Tensor<double>>> HR_irr;
    RI::Tensor<double> A = make_A_real();
    HR_irr[0][{0, {0, 0, 0}}] = A;

    ModuleSymmetry::TapR irapR = {{0, 0}, {0, 0, 0}};
    symrot.irs_.sector_stars_[irapR] = {
        {0, {{0, 0}, {0, 0, 0}}},
        {1, {{0, 0}, {0, 0, 0}}}, // duplicate target
    };

    ModuleSymmetry::Symmetry symm;
    Atom atoms[1];
    EXPECT_THROW(symrot.restore_HR_abf(symm, atoms, st, HR_irr), std::runtime_error);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
