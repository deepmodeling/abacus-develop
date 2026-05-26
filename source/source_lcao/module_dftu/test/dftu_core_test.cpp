#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include <numeric>
#include <algorithm>

/***********************************************************************
 * Unit tests for DFT+U core algorithms.
 *
 * These tests target the most complex and bug-prone logic:
 * 1. eff_pot_pw_index calculation for mixed atom types and nspin modes
 * 2. copy_locale <-> set_locale roundtrip (3 data layouts)
 * 3. VU effective potential formula (cal_type=3, FLL)
 * 4. nspin=1 symmetrization: (n + n^T) / 2
 * 5. Locale mixing matrix operations
 ***********************************************************************/

// =====================================================================
// 1. eff_pot_pw_index calculation
//
// nspin=1: offset = sum(tlp1^2), total = sum(all tlp1^2)
// nspin=2: same per-spin-channel, then pot_index *= 2 (split layout)
// nspin=4: offset = sum((tlp1*npol)^2), each atom = 4*tlp1^2
// =====================================================================

class EffPotIndexTest : public ::testing::Test
{
  protected:
    struct AtomSpec { int l; int na; }; // correlated orbital l, number of atoms
    std::vector<int> eff_pot_pw_index;
    int pot_index;

    void compute_indices(const std::vector<AtomSpec>& atoms, int nspin)
    {
        pot_index = 0;
        eff_pot_pw_index.resize(atoms.size());

        for (size_t i = 0; i < atoms.size(); i++)
        {
            int tlp1 = 2 * atoms[i].l + 1;
            int tlp1_npol = tlp1 * (nspin == 4 ? 2 : 1);

            if (nspin == 4)
            {
                eff_pot_pw_index[i] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                eff_pot_pw_index[i] = pot_index;
                pot_index += tlp1 * tlp1;
            }
        }

        if (nspin == 2)
            pot_index *= 2;
    }
};

TEST_F(EffPotIndexTest, Nspin1_MixedOrbitals)
{
    // 3 atoms: p(l=1), d(l=2), p(l=1)
    std::vector<AtomSpec> atoms = {{1, 1}, {2, 1}, {1, 1}};
    compute_indices(atoms, 1);

    // p: 9, d: 25, p: 9
    EXPECT_EQ(eff_pot_pw_index[0], 0);
    EXPECT_EQ(eff_pot_pw_index[1], 9);
    EXPECT_EQ(eff_pot_pw_index[2], 34);
    EXPECT_EQ(pot_index, 43); // 9 + 25 + 9
}

TEST_F(EffPotIndexTest, Nspin2_SplitLayout)
{
    // 2 atoms: d(l=2), d(l=2)
    std::vector<AtomSpec> atoms = {{2, 1}, {2, 1}};
    compute_indices(atoms, 2);

    // Each atom: 25 entries per spin channel
    // Split: [atom0_up, atom1_up | atom0_dn, atom1_dn]
    EXPECT_EQ(eff_pot_pw_index[0], 0);
    EXPECT_EQ(eff_pot_pw_index[1], 25);
    EXPECT_EQ(pot_index, 100); // (25 + 25) * 2
}

TEST_F(EffPotIndexTest, Nspin4_PauliBlocks)
{
    // 2 atoms: d(l=2), p(l=1)
    std::vector<AtomSpec> atoms = {{2, 1}, {1, 1}};
    compute_indices(atoms, 4);

    // nspin=4: tlp1_npol = tlp1 * 2, size = (tlp1*2)^2 = 4*tlp1^2
    // atom0 (d): (5*2)^2 = 100
    // atom1 (p): (3*2)^2 = 36
    EXPECT_EQ(eff_pot_pw_index[0], 0);
    EXPECT_EQ(eff_pot_pw_index[1], 100);
    EXPECT_EQ(pot_index, 136); // 100 + 36
}

TEST_F(EffPotIndexTest, Nspin2_MultiAtomSameType_Ordering)
{
    // 4 atoms: all d(l=2), testing ordering consistency
    std::vector<AtomSpec> atoms = {{2, 1}, {2, 1}, {2, 1}, {2, 1}};
    compute_indices(atoms, 2);

    for (int i = 0; i < 4; i++)
        EXPECT_EQ(eff_pot_pw_index[i], i * 25);

    EXPECT_EQ(pot_index, 200); // 4 * 25 * 2
}

TEST_F(EffPotIndexTest, Nspin1_LargeSystem)
{
    // 10 atoms: mixed p and d
    std::vector<AtomSpec> atoms;
    for (int i = 0; i < 10; i++)
        atoms.push_back({(i % 2 == 0) ? 2 : 1, 1});

    compute_indices(atoms, 1);

    int expected = 0;
    for (int i = 0; i < 10; i++)
    {
        EXPECT_EQ(eff_pot_pw_index[i], expected);
        int tlp1 = 2 * atoms[i].l + 1;
        expected += tlp1 * tlp1;
    }
    EXPECT_EQ(pot_index, expected);
}

// =====================================================================
// 2. copy_locale <-> set_locale roundtrip
//
// Tests the bidirectional conversion between nested locale matrix
// and flat uom_array/uom_save arrays for all 3 nspin modes.
// =====================================================================

struct Matrix2D {
    int nr, nc;
    std::vector<double> data;
    Matrix2D() : nr(0), nc(0), data() {}
    Matrix2D(int r, int c) : nr(r), nc(c), data(r * c, 0.0) {}
    double& operator()(int i, int j) { return data[i * nc + j]; }
    const double& operator()(int i, int j) const { return data[i * nc + j]; }
};

static void copy_locale_to_flat(
    const std::vector<Matrix2D>& locale_up,
    const std::vector<Matrix2D>& locale_dn,
    std::vector<double>& uom_save,
    const std::vector<int>& eff_pot_pw_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[eff_pot_pw_index[iat] + mm] = locale_up[iat].data[mm];
        }
    }
    else if (nspin == 2) // split layout: [up | dn]
    {
        int half_size = uom_save.size() / 2;
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                uom_save[eff_pot_pw_index[iat] + mm] = locale_up[iat].data[mm];
                uom_save[half_size + eff_pot_pw_index[iat] + mm] = locale_dn[iat].data[mm];
            }
        }
    }
    else // nspin=1: single spin channel
    {
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[eff_pot_pw_index[iat] + mm] = locale_up[iat].data[mm];
        }
    }
}

static void set_locale_from_flat(
    const std::vector<double>& uom_array,
    std::vector<Matrix2D>& locale_up,
    std::vector<Matrix2D>& locale_dn,
    const std::vector<int>& eff_pot_pw_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                locale_up[iat].data[mm] = uom_array[eff_pot_pw_index[iat] + mm];
        }
    }
    else if (nspin == 2)
    {
        int half_size = uom_array.size() / 2;
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                locale_up[iat].data[mm] = uom_array[eff_pot_pw_index[iat] + mm];
                locale_dn[iat].data[mm] = uom_array[half_size + eff_pot_pw_index[iat] + mm];
            }
        }
    }
    else // nspin=1
    {
        for (size_t iat = 0; iat < locale_up.size(); iat++)
        {
            int size = locale_up[iat].nr * locale_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                locale_up[iat].data[mm] = uom_array[eff_pot_pw_index[iat] + mm];
        }
    }
}

class LocaleRoundtripTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(LocaleRoundtripTest, Nspin1_SingleAtom)
{
    const int l = 2; // d-orbital
    const int size = (2 * l + 1) * (2 * l + 1); // 25
    int nspin = 1;

    std::vector<Matrix2D> locale_up(1, Matrix2D(2 * l + 1, 2 * l + 1));
    std::vector<Matrix2D> locale_dn(1, Matrix2D(2 * l + 1, 2 * l + 1));
    for (int i = 0; i < size; i++)
    {
        locale_up[0].data[i] = static_cast<double>(i + 1);
        locale_dn[0].data[i] = static_cast<double>(i + 1); // nspin=1: same
    }

    std::vector<int> eff_pot_pw_index = {0};
    std::vector<double> uom_save(size, 0.0);

    copy_locale_to_flat(locale_up, locale_dn, uom_save, eff_pot_pw_index, nspin);
    set_locale_from_flat(uom_save, locale_up, locale_dn, eff_pot_pw_index, nspin);

    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(locale_up[0].data[i], static_cast<double>(i + 1));
}

TEST_F(LocaleRoundtripTest, Nspin2_SplitLayout)
{
    const int l = 2;
    const int size = (2 * l + 1) * (2 * l + 1); // 25
    const int total = size * 2; // split: [up | dn]
    int nspin = 2;

    std::vector<Matrix2D> locale_up(1, Matrix2D(2 * l + 1, 2 * l + 1));
    std::vector<Matrix2D> locale_dn(1, Matrix2D(2 * l + 1, 2 * l + 1));
    for (int i = 0; i < size; i++)
    {
        locale_up[0].data[i] = static_cast<double>(i + 1);
        locale_dn[0].data[i] = static_cast<double>(i + 100);
    }

    std::vector<int> eff_pot_pw_index = {0};
    std::vector<double> uom_save(total, 0.0);

    copy_locale_to_flat(locale_up, locale_dn, uom_save, eff_pot_pw_index, nspin);

    // Verify split layout: first half = up, second half = dn
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(uom_save[size + i], static_cast<double>(i + 100));
    }

    set_locale_from_flat(uom_save, locale_up, locale_dn, eff_pot_pw_index, nspin);

    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(locale_up[0].data[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(locale_dn[0].data[i], static_cast<double>(i + 100));
    }
}

TEST_F(LocaleRoundtripTest, Nspin2_MultiAtom_SplitLayout)
{
    // 3 atoms with different orbital sizes
    struct AtomSpec { int l; int start_val; };
    std::vector<AtomSpec> specs = {{1, 1}, {2, 100}, {1, 1000}}; // p, d, p
    int nspin = 2;

    std::vector<int> sizes;
    for (auto& s : specs)
        sizes.push_back((2 * s.l + 1) * (2 * s.l + 1));

    int half_size = std::accumulate(sizes.begin(), sizes.end(), 0);
    int total = half_size * 2;

    std::vector<int> eff_pot_pw_index(specs.size());
    int offset = 0;
    for (size_t i = 0; i < specs.size(); i++)
    {
        eff_pot_pw_index[i] = offset;
        offset += sizes[i];
    }

    std::vector<Matrix2D> locale_up(specs.size()), locale_dn(specs.size());
    for (size_t i = 0; i < specs.size(); i++)
    {
        int dim = 2 * specs[i].l + 1;
        locale_up[i] = Matrix2D(dim, dim);
        locale_dn[i] = Matrix2D(dim, dim);
        for (int j = 0; j < sizes[i]; j++)
        {
            locale_up[i].data[j] = static_cast<double>(specs[i].start_val + j);
            locale_dn[i].data[j] = static_cast<double>(specs[i].start_val + j + 5000);
        }
    }

    std::vector<double> uom_save(total, 0.0);
    copy_locale_to_flat(locale_up, locale_dn, uom_save, eff_pot_pw_index, nspin);
    set_locale_from_flat(uom_save, locale_up, locale_dn, eff_pot_pw_index, nspin);

    for (size_t i = 0; i < specs.size(); i++)
    {
        for (int j = 0; j < sizes[i]; j++)
        {
            EXPECT_DOUBLE_EQ(locale_up[i].data[j], static_cast<double>(specs[i].start_val + j));
            EXPECT_DOUBLE_EQ(locale_dn[i].data[j], static_cast<double>(specs[i].start_val + j + 5000));
        }
    }
}

TEST_F(LocaleRoundtripTest, Nspin4_PauliBlocks)
{
    // 2 atoms: d(l=2), p(l=1)
    struct AtomSpec { int l; };
    std::vector<AtomSpec> specs = {{2}, {1}};
    int nspin = 4;
    int npol = 2;

    std::vector<int> sizes;
    for (auto& s : specs)
    {
        int tlp1 = 2 * s.l + 1;
        sizes.push_back((tlp1 * npol) * (tlp1 * npol)); // 4*tlp1^2
    }

    int total = std::accumulate(sizes.begin(), sizes.end(), 0);

    std::vector<int> eff_pot_pw_index(specs.size());
    int offset = 0;
    for (size_t i = 0; i < specs.size(); i++)
    {
        eff_pot_pw_index[i] = offset;
        offset += sizes[i];
    }

    std::vector<Matrix2D> locale(specs.size());
    for (size_t i = 0; i < specs.size(); i++)
    {
        int dim = (2 * specs[i].l + 1) * npol;
        locale[i] = Matrix2D(dim, dim);
        for (int j = 0; j < sizes[i]; j++)
            locale[i].data[j] = static_cast<double>(i * 1000 + j + 1);
    }

    std::vector<double> uom_array(total, 0.0);
    std::vector<Matrix2D> locale_dn(specs.size()); // unused for nspin=4

    // nspin=4: only locale_up is used
    copy_locale_to_flat(locale, locale_dn, uom_array, eff_pot_pw_index, nspin);
    set_locale_from_flat(uom_array, locale, locale_dn, eff_pot_pw_index, nspin);

    for (size_t i = 0; i < specs.size(); i++)
    {
        for (int j = 0; j < sizes[i]; j++)
            EXPECT_DOUBLE_EQ(locale[i].data[j], static_cast<double>(i * 1000 + j + 1));
    }
}

// =====================================================================
// 3. VU effective potential formula (cal_type=3, FLL)
//
// VU[m0,m1] = U * (0.5*delta(m0,m1) - locale[m0,m1])  (diagonal)
// VU[m0,m1] = -U * locale[m0,m1]                       (off-diagonal)
// =====================================================================

static double compute_vu(double U_val, int m0, int m1, double locale_val)
{
    if (m0 == m1)
        return U_val * (0.5 - locale_val);
    else
        return -U_val * locale_val;
}

class VUPotentialTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(VUPotentialTest, Diagonal_HalfFilled)
{
    double U = 4.0;
    double locale = 0.5; // half-filled
    double vu = compute_vu(U, 0, 0, locale);
    EXPECT_DOUBLE_EQ(vu, 0.0); // U * (0.5 - 0.5) = 0
}

TEST_F(VUPotentialTest, Diagonal_FullyOccupied)
{
    double U = 4.0;
    double locale = 1.0; // fully occupied
    double vu = compute_vu(U, 0, 0, locale);
    EXPECT_DOUBLE_EQ(vu, -2.0); // U * (0.5 - 1.0) = -2.0
}

TEST_F(VUPotentialTest, Diagonal_Empty)
{
    double U = 4.0;
    double locale = 0.0; // empty
    double vu = compute_vu(U, 0, 0, locale);
    EXPECT_DOUBLE_EQ(vu, 2.0); // U * 0.5 = 2.0
}

TEST_F(VUPotentialTest, OffDiagonal)
{
    double U = 5.0;
    double locale = 0.3;
    double vu = compute_vu(U, 0, 1, locale);
    EXPECT_DOUBLE_EQ(vu, -1.5); // -U * locale = -1.5
}

TEST_F(VUPotentialTest, FullMatrix_DOrbital)
{
    double U = 4.0;
    const int m_size = 5; // d-orbital

    // Create a locale matrix with known values
    std::vector<double> locale(m_size * m_size, 0.0);
    for (int m = 0; m < m_size; m++)
        locale[m * m_size + m] = 0.8 - 0.1 * m; // diagonal varies

    // Compute VU matrix
    std::vector<double> vu(m_size * m_size, 0.0);
    for (int m0 = 0; m0 < m_size; m0++)
        for (int m1 = 0; m1 < m_size; m1++)
            vu[m0 * m_size + m1] = compute_vu(U, m0, m1, locale[m0 * m_size + m1]);

    // Verify diagonal
    for (int m = 0; m < m_size; m++)
    {
        double expected = U * (0.5 - (0.8 - 0.1 * m));
        EXPECT_DOUBLE_EQ(vu[m * m_size + m], expected);
    }

    // Verify off-diagonal (should be zero since locale is diagonal)
    for (int m0 = 0; m0 < m_size; m0++)
        for (int m1 = 0; m1 < m_size; m1++)
            if (m0 != m1)
                EXPECT_DOUBLE_EQ(vu[m0 * m_size + m1], 0.0);
}

TEST_F(VUPotentialTest, WithYukawa_UminusJ)
{
    // With Yukawa: use (U - J) instead of U
    double U_val = 5.0, J_val = 1.0;
    double U_eff = U_val - J_val;

    double vu_diag = compute_vu(U_eff, 0, 0, 0.6);
    double vu_off = compute_vu(U_eff, 0, 1, 0.2);

    EXPECT_DOUBLE_EQ(vu_diag, U_eff * (0.5 - 0.6)); // -0.4
    EXPECT_DOUBLE_EQ(vu_off, -U_eff * 0.2); // -0.8
}

// =====================================================================
// 4. nspin=1 symmetrization: (n + n^T) / 2
//
// After MPI_Allreduce, nspin=1 locale is symmetrized:
//   locale = (locale + locale^T) / 2
//   locale[1] = locale[0]  (copy to spin-down)
// =====================================================================

static void symmetrize_nspin1(std::vector<double>& locale, int m_size)
{
    std::vector<double> temp = locale;
    for (int m0 = 0; m0 < m_size; m0++)
        for (int m1 = 0; m1 < m_size; m1++)
            locale[m0 * m_size + m1] = (temp[m0 * m_size + m1] + temp[m1 * m_size + m0]) / 2.0;
}

class NSpin1SymmetrizationTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(NSpin1SymmetrizationTest, AlreadySymmetric)
{
    const int m_size = 3;
    std::vector<double> locale = {
        0.5, 0.1, 0.0,
        0.1, 0.3, 0.2,
        0.0, 0.2, 0.2
    };
    std::vector<double> expected = locale;

    symmetrize_nspin1(locale, m_size);

    for (int i = 0; i < m_size * m_size; i++)
        EXPECT_DOUBLE_EQ(locale[i], expected[i]);
}

TEST_F(NSpin1SymmetrizationTest, AsymmetricInput)
{
    const int m_size = 2;
    std::vector<double> locale = {
        0.5, 0.3,
        0.1, 0.5
    };

    symmetrize_nspin1(locale, m_size);

    // locale[0,1] = (0.3 + 0.1) / 2 = 0.2
    // locale[1,0] = (0.1 + 0.3) / 2 = 0.2
    EXPECT_DOUBLE_EQ(locale[0], 0.5);
    EXPECT_DOUBLE_EQ(locale[1], 0.2);
    EXPECT_DOUBLE_EQ(locale[2], 0.2);
    EXPECT_DOUBLE_EQ(locale[3], 0.5);
}

TEST_F(NSpin1SymmetrizationTest, ResultIsSymmetric)
{
    const int m_size = 5;
    std::vector<double> locale(m_size * m_size);
    for (int i = 0; i < m_size * m_size; i++)
        locale[i] = static_cast<double>(rand()) / RAND_MAX; // random values

    symmetrize_nspin1(locale, m_size);

    // Verify symmetry
    for (int m0 = 0; m0 < m_size; m0++)
        for (int m1 = 0; m1 < m_size; m1++)
            EXPECT_DOUBLE_EQ(locale[m0 * m_size + m1], locale[m1 * m_size + m0]);
}

TEST_F(NSpin1SymmetrizationTest, TracePreserved)
{
    const int m_size = 3;
    std::vector<double> locale = {
        0.5, 0.1, 0.2,
        0.3, 0.4, 0.1,
        0.2, 0.1, 0.3
    };

    double trace_before = locale[0] + locale[4] + locale[8];
    symmetrize_nspin1(locale, m_size);
    double trace_after = locale[0] + locale[4] + locale[8];

    EXPECT_DOUBLE_EQ(trace_before, trace_after);
}

// =====================================================================
// 5. Locale mixing matrix operations
//
// locale_new = locale_current * beta + locale_save * (1 - beta)
// =====================================================================

class LocaleMixingTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(LocaleMixingTest, MixingWithBeta)
{
    const int size = 25;
    std::vector<double> locale_current(size), locale_save(size);
    for (int i = 0; i < size; i++)
    {
        locale_current[i] = static_cast<double>(i + 1);
        locale_save[i] = static_cast<double>(i + 100);
    }

    double beta = 0.3;
    std::vector<double> locale_new(size);
    for (int i = 0; i < size; i++)
        locale_new[i] = locale_current[i] * beta + locale_save[i] * (1.0 - beta);

    // Spot check
    EXPECT_DOUBLE_EQ(locale_new[0], 1.0 * 0.3 + 100.0 * 0.7);
    EXPECT_DOUBLE_EQ(locale_new[12], 13.0 * 0.3 + 112.0 * 0.7);
    EXPECT_DOUBLE_EQ(locale_new[24], 25.0 * 0.3 + 124.0 * 0.7);
}

TEST_F(LocaleMixingTest, BetaZero_NoChange)
{
    const int size = 9;
    std::vector<double> locale_current(size), locale_save(size);
    for (int i = 0; i < size; i++)
    {
        locale_current[i] = static_cast<double>(i);
        locale_save[i] = static_cast<double>(i + 50);
    }

    double beta = 0.0;
    std::vector<double> locale_new(size);
    for (int i = 0; i < size; i++)
        locale_new[i] = locale_current[i] * beta + locale_save[i] * (1.0 - beta);

    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(locale_new[i], locale_save[i]);
}

TEST_F(LocaleMixingTest, BetaOne_FullReplace)
{
    const int size = 9;
    std::vector<double> locale_current(size), locale_save(size);
    for (int i = 0; i < size; i++)
    {
        locale_current[i] = static_cast<double>(i);
        locale_save[i] = static_cast<double>(i + 50);
    }

    double beta = 1.0;
    std::vector<double> locale_new(size);
    for (int i = 0; i < size; i++)
        locale_new[i] = locale_current[i] * beta + locale_save[i] * (1.0 - beta);

    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(locale_new[i], locale_current[i]);
}

TEST_F(LocaleMixingTest, MatrixMixing_Nspin2)
{
    // nspin=2: mix both spin channels independently
    const int size = 25;
    std::vector<double> current_up(size), current_dn(size);
    std::vector<double> save_up(size), save_dn(size);

    for (int i = 0; i < size; i++)
    {
        current_up[i] = static_cast<double>(i);
        current_dn[i] = static_cast<double>(i + 1000);
        save_up[i] = static_cast<double>(i + 100);
        save_dn[i] = static_cast<double>(i + 1100);
    }

    double beta = 0.4;
    std::vector<double> new_up(size), new_dn(size);
    for (int i = 0; i < size; i++)
    {
        new_up[i] = current_up[i] * beta + save_up[i] * (1.0 - beta);
        new_dn[i] = current_dn[i] * beta + save_dn[i] * (1.0 - beta);
    }

    EXPECT_DOUBLE_EQ(new_up[0], 0.0 * 0.4 + 100.0 * 0.6);
    EXPECT_DOUBLE_EQ(new_dn[0], 1000.0 * 0.4 + 1100.0 * 0.6);
}

// =====================================================================
// 6. Uramping logic
// =====================================================================

class UrampingTest : public ::testing::Test
{
  protected:
    std::vector<double> U;
    std::vector<double> U0;
    double uramping;

    void SetUp() override
    {
        U = {0.0, 0.0, 0.0};
        U0 = {4.0, 6.0, 8.0};
        uramping = 1.0;
    }

    void apply_uramping()
    {
        if (uramping < 0.01) return;
        for (size_t i = 0; i < U0.size(); i++)
        {
            if (U[i] + uramping < U0[i])
                U[i] += uramping;
            else
                U[i] = U0[i];
        }
    }

    bool u_converged()
    {
        for (size_t i = 0; i < U0.size(); i++)
            if (U[i] != U0[i]) return false;
        return true;
    }
};

TEST_F(UrampingTest, RampUp_DifferentRates)
{
    for (int step = 0; step < 10; step++)
        apply_uramping();

    EXPECT_DOUBLE_EQ(U[0], 4.0); // capped at step 4
    EXPECT_DOUBLE_EQ(U[1], 6.0); // capped at step 6
    EXPECT_DOUBLE_EQ(U[2], 8.0); // capped at step 8
    EXPECT_TRUE(u_converged());
}

TEST_F(UrampingTest, UrampingDisabled)
{
    uramping = -1.0;
    apply_uramping();
    EXPECT_DOUBLE_EQ(U[0], 0.0);
    EXPECT_DOUBLE_EQ(U[1], 0.0);
    EXPECT_DOUBLE_EQ(U[2], 0.0);
}

TEST_F(UrampingTest, PrecisionIssue_DoubleComparison)
{
    // Test potential floating-point comparison issues
    U = {0.0};
    U0 = {1.0 / 3.0}; // 0.333...
    uramping = 1.0 / 3.0;

    apply_uramping(); // Should reach U0
    // Note: exact equality may fail due to floating-point
    EXPECT_NEAR(U[0], U0[0], 1e-15);
}

// =====================================================================
// 7. Energy correction formula
//
// E_U = 0.5 * U * sum_spin [Tr(n) - Tr(n^2)]
// =====================================================================

class EnergyCorrectionTest : public ::testing::Test
{
  protected:
    static double compute_energy(const std::vector<double>& locale_flat, int m_size, double U)
    {
        double nm_trace = 0.0, nm2_trace = 0.0;
        for (int m0 = 0; m0 < m_size; m0++)
        {
            nm_trace += locale_flat[m0 * m_size + m0];
            for (int m1 = 0; m1 < m_size; m1++)
                nm2_trace += locale_flat[m0 * m_size + m1] * locale_flat[m1 * m_size + m0];
        }
        return 0.5 * U * (nm_trace - nm2_trace);
    }
};

TEST_F(EnergyCorrectionTest, HalfFilled_DOrbital)
{
    const int m_size = 5;
    std::vector<double> locale(m_size * m_size, 0.0);
    for (int m = 0; m < m_size; m++)
        locale[m * m_size + m] = 0.5;

    double energy = compute_energy(locale, m_size, 4.0);
    // Tr(n) = 2.5, Tr(n^2) = 1.25, E = 0.5 * 4 * 1.25 = 2.5
    EXPECT_DOUBLE_EQ(energy, 2.5);
}

TEST_F(EnergyCorrectionTest, FullyOccupied_ZeroEnergy)
{
    const int m_size = 5;
    std::vector<double> locale(m_size * m_size, 0.0);
    for (int m = 0; m < m_size; m++)
        locale[m * m_size + m] = 1.0;

    double energy = compute_energy(locale, m_size, 4.0);
    // Tr(n) = Tr(n^2) = 5, E = 0
    EXPECT_DOUBLE_EQ(energy, 0.0);
}

TEST_F(EnergyCorrectionTest, Empty_ZeroEnergy)
{
    const int m_size = 5;
    std::vector<double> locale(m_size * m_size, 0.0);
    double energy = compute_energy(locale, m_size, 4.0);
    EXPECT_DOUBLE_EQ(energy, 0.0);
}

TEST_F(EnergyCorrectionTest, OffDiagonal_Contribution)
{
    const int m_size = 2;
    std::vector<double> locale = {
        0.3, 0.1,
        0.1, 0.3
    };

    double energy = compute_energy(locale, m_size, 4.0);
    // Tr(n) = 0.6, Tr(n^2) = 0.3^2 + 0.1^2 + 0.1^2 + 0.3^2 = 0.20
    // E = 0.5 * 4 * (0.6 - 0.20) = 0.8
    EXPECT_DOUBLE_EQ(energy, 0.8);
}

TEST_F(EnergyCorrectionTest, DoubleCounting_Energy)
{
    // E_dc = sum_{m1,m2,spin} VU[m1,m2] * n[m2,m1]
    const int m_size = 3;
    double U = 4.0;
    std::vector<double> locale = {
        0.5, 0.0, 0.0,
        0.0, 0.3, 0.0,
        0.0, 0.0, 0.2
    };

    double e_dc = 0.0;
    for (int m1 = 0; m1 < m_size; m1++)
        for (int m2 = 0; m2 < m_size; m2++)
        {
            double vu = (m1 == m2) ? U * (0.5 - locale[m1 * m_size + m2])
                                   : -U * locale[m1 * m_size + m2];
            e_dc += vu * locale[m2 * m_size + m1];
        }

    // Only diagonal: m=0: 0*0.5=0, m=1: 0.8*0.3=0.24, m=2: 1.2*0.2=0.24
    EXPECT_NEAR(e_dc, 0.48, 1e-14);
}

// =====================================================================
// 8. Memory layout and boundary conditions
// =====================================================================

class MemoryLayoutTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(MemoryLayoutTest, FlatArray_ContiguousMemory)
{
    // Verify that flat arrays are contiguous for efficient memory access
    const int n_atoms = 100;
    const int l = 2;
    const int size = (2 * l + 1) * (2 * l + 1);

    std::vector<double> uom_array(n_atoms * size * 2, 0.0); // nspin=2

    // Fill with stride-1 access pattern (cache-friendly)
    for (size_t i = 0; i < uom_array.size(); i++)
        uom_array[i] = static_cast<double>(i);

    // Verify sequential access works correctly
    double sum = 0.0;
    for (size_t i = 0; i < uom_array.size(); i++)
        sum += uom_array[i];

    double expected = uom_array.size() * (uom_array.size() - 1) / 2.0;
    EXPECT_DOUBLE_EQ(sum, expected);
}

TEST_F(MemoryLayoutTest, IndexBounds_MultiAtom)
{
    // Test that indices don't exceed array bounds
    std::vector<int> atom_sizes = {9, 25, 9, 25, 9}; // p, d, p, d, p
    int total = std::accumulate(atom_sizes.begin(), atom_sizes.end(), 0);

    std::vector<int> eff_pot_pw_index(atom_sizes.size());
    int offset = 0;
    for (size_t i = 0; i < atom_sizes.size(); i++)
    {
        eff_pot_pw_index[i] = offset;
        offset += atom_sizes[i];
    }

    // nspin=2: double the size
    int total_split = total * 2;
    std::vector<double> uom_array(total_split, 0.0);

    // Verify all indices are within bounds
    for (size_t i = 0; i < atom_sizes.size(); i++)
    {
        int start = eff_pot_pw_index[i];
        int end = start + atom_sizes[i];
        EXPECT_GE(start, 0);
        EXPECT_LE(end, total); // within half_size

        // Also check the spin-down half
        EXPECT_GE(total + start, total);
        EXPECT_LE(total + end, total_split);
    }
}

TEST_F(MemoryLayoutTest, LargeMatrix_NumericalStability)
{
    // Test numerical stability with large matrices
    const int m_size = 7; // f-orbital
    const int size = m_size * m_size;

    std::vector<double> locale(size, 0.0);
    // Create a symmetric matrix with diagonal elements in [0,1]
    // This represents a valid occupation matrix
    for (int m = 0; m < m_size; m++)
        locale[m * m_size + m] = 0.5 + (static_cast<double>(rand()) / RAND_MAX - 0.5) * 0.001;

    // Symmetrize: n = (n + n^T) / 2
    std::vector<double> symm_locale(size, 0.0);
    for (int m0 = 0; m0 < m_size; m0++)
        for (int m1 = 0; m1 < m_size; m1++)
            symm_locale[m0 * m_size + m1] = (locale[m0 * m_size + m1] + locale[m1 * m_size + m0]) / 2.0;

    double nm_trace = 0.0, nm2_trace = 0.0;
    for (int m0 = 0; m0 < m_size; m0++)
    {
        nm_trace += symm_locale[m0 * m_size + m0];
        for (int m1 = 0; m1 < m_size; m1++)
            nm2_trace += symm_locale[m0 * m_size + m1] * symm_locale[m1 * m_size + m0];
    }

    // For half-filled symmetric matrix, Tr(n) - Tr(n^2) should be positive
    double diff = nm_trace - nm2_trace;
    EXPECT_GT(diff, 0.0);
    EXPECT_LT(diff, static_cast<double>(m_size)); // bounded by m_size
}
