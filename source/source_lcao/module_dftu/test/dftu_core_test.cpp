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
 * 4. Energy correction and double-counting terms
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

TEST_F(EffPotIndexTest, Nspin2and4_SplitAndPauli)
{
    // nspin=2: 2 d-atoms, split layout [up | dn]
    std::vector<AtomSpec> atoms2 = {{2, 1}, {2, 1}};
    compute_indices(atoms2, 2);
    EXPECT_EQ(eff_pot_pw_index[0], 0);
    EXPECT_EQ(eff_pot_pw_index[1], 25);
    EXPECT_EQ(pot_index, 100); // (25 + 25) * 2

    // nspin=4: d + p atoms, Pauli blocks
    std::vector<AtomSpec> atoms4 = {{2, 1}, {1, 1}};
    compute_indices(atoms4, 4);
    EXPECT_EQ(eff_pot_pw_index[0], 0);    // d: (5*2)^2 = 100
    EXPECT_EQ(eff_pot_pw_index[1], 100);  // p: (3*2)^2 = 36
    EXPECT_EQ(pot_index, 136);
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

TEST_F(LocaleRoundtripTest, Nspin1and2_SingleAndSplitLayout)
{
    // nspin=1: single atom d-orbital roundtrip
    const int l = 2;
    const int size = (2 * l + 1) * (2 * l + 1); // 25

    std::vector<Matrix2D> locale_up(1, Matrix2D(2 * l + 1, 2 * l + 1));
    std::vector<Matrix2D> locale_dn(1, Matrix2D(2 * l + 1, 2 * l + 1));
    for (int i = 0; i < size; i++)
        locale_up[0].data[i] = static_cast<double>(i + 1);

    std::vector<int> eff_pot_pw_index = {0};
    std::vector<double> uom_save(size, 0.0);
    copy_locale_to_flat(locale_up, locale_dn, uom_save, eff_pot_pw_index, 1);
    set_locale_from_flat(uom_save, locale_up, locale_dn, eff_pot_pw_index, 1);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(locale_up[0].data[i], static_cast<double>(i + 1));

    // nspin=2: split layout [up | dn] with distinct values
    const int total = size * 2;
    for (int i = 0; i < size; i++)
    {
        locale_up[0].data[i] = static_cast<double>(i + 1);
        locale_dn[0].data[i] = static_cast<double>(i + 100);
    }
    uom_save.assign(total, 0.0);
    copy_locale_to_flat(locale_up, locale_dn, uom_save, eff_pot_pw_index, 2);
    // Verify split layout
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(uom_save[size + i], static_cast<double>(i + 100));
    }
    set_locale_from_flat(uom_save, locale_up, locale_dn, eff_pot_pw_index, 2);
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(locale_up[0].data[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(locale_dn[0].data[i], static_cast<double>(i + 100));
    }
}

TEST_F(LocaleRoundtripTest, Nspin4_PauliBlocks)
{
    // 2 atoms: d(l=2), p(l=1)
    struct AtomSpec { int l; };
    std::vector<AtomSpec> specs = {{2}, {1}};
    int npol = 2;

    std::vector<int> sizes;
    for (auto& s : specs)
    {
        int tlp1 = 2 * s.l + 1;
        sizes.push_back((tlp1 * npol) * (tlp1 * npol));
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

    copy_locale_to_flat(locale, locale_dn, uom_array, eff_pot_pw_index, 4);
    set_locale_from_flat(uom_array, locale, locale_dn, eff_pot_pw_index, 4);

    for (size_t i = 0; i < specs.size(); i++)
        for (int j = 0; j < sizes[i]; j++)
            EXPECT_DOUBLE_EQ(locale[i].data[j], static_cast<double>(i * 1000 + j + 1));
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

TEST_F(VUPotentialTest, OffDiagonal)
{
    double U = 5.0;
    double locale = 0.3;
    double vu = compute_vu(U, 0, 1, locale);
    EXPECT_DOUBLE_EQ(vu, -1.5); // -U * locale = -1.5
}

// =====================================================================
// 4. Energy correction formula
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
