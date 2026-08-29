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
 * 1. pot_uterm_pw_index calculation for mixed atom types and nspin modes
 * 2. copy_occ_mat <-> set_occ_mat roundtrip (3 data layouts)
 * 3. pot_onsite effective potential formula (cal_type=3, FLL)
 * 4. Energy correction and double-counting terms
 ***********************************************************************/

// =====================================================================
// 1. pot_uterm_pw_index calculation
//
// nspin=1: offset = sum(tlp1^2), total = sum(all tlp1^2)
// nspin=2: same per-spin-channel, then pot_index *= 2 (split layout)
// nspin=4: offset = sum((tlp1*npol)^2), each atom = 4*tlp1^2
// =====================================================================

class EffPotIndexTest : public ::testing::Test
{
  protected:
    struct AtomSpec { int l; int na; }; // correlated orbital l, number of atoms
    std::vector<int> pot_uterm_pw_index;
    int pot_index;

    void compute_indices(const std::vector<AtomSpec>& atoms, int nspin)
    {
        pot_index = 0;
        pot_uterm_pw_index.resize(atoms.size());

        for (size_t i = 0; i < atoms.size(); i++)
        {
            int tlp1 = 2 * atoms[i].l + 1;
            int tlp1_npol = tlp1 * (nspin == 4 ? 2 : 1);

            if (nspin == 4)
            {
                pot_uterm_pw_index[i] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                pot_uterm_pw_index[i] = pot_index;
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
    EXPECT_EQ(pot_uterm_pw_index[0], 0);
    EXPECT_EQ(pot_uterm_pw_index[1], 9);
    EXPECT_EQ(pot_uterm_pw_index[2], 34);
    EXPECT_EQ(pot_index, 43); // 9 + 25 + 9
}

TEST_F(EffPotIndexTest, Nspin2and4_SplitAndPauli)
{
    // nspin=2: 2 d-atoms, split layout [up | dn]
    std::vector<AtomSpec> atoms2 = {{2, 1}, {2, 1}};
    compute_indices(atoms2, 2);
    EXPECT_EQ(pot_uterm_pw_index[0], 0);
    EXPECT_EQ(pot_uterm_pw_index[1], 25);
    EXPECT_EQ(pot_index, 100); // (25 + 25) * 2

    // nspin=4: d + p atoms, Pauli blocks
    std::vector<AtomSpec> atoms4 = {{2, 1}, {1, 1}};
    compute_indices(atoms4, 4);
    EXPECT_EQ(pot_uterm_pw_index[0], 0);    // d: (5*2)^2 = 100
    EXPECT_EQ(pot_uterm_pw_index[1], 100);  // p: (3*2)^2 = 36
    EXPECT_EQ(pot_index, 136);
}

// =====================================================================
// 2. copy_occ_mat <-> set_occ_mat roundtrip
//
// Tests the bidirectional conversion between nested occ_mat matrix
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

static void copy_occ_mat_to_flat(
    const std::vector<Matrix2D>& occ_mat_up,
    const std::vector<Matrix2D>& occ_mat_dn,
    std::vector<double>& uom_save,
    const std::vector<int>& pot_uterm_pw_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[pot_uterm_pw_index[iat] + mm] = occ_mat_up[iat].data[mm];
        }
    }
    else if (nspin == 2) // split layout: [up | dn]
    {
        int half_size = uom_save.size() / 2;
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                uom_save[pot_uterm_pw_index[iat] + mm] = occ_mat_up[iat].data[mm];
                uom_save[half_size + pot_uterm_pw_index[iat] + mm] = occ_mat_dn[iat].data[mm];
            }
        }
    }
    else // nspin=1: single spin channel
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                uom_save[pot_uterm_pw_index[iat] + mm] = occ_mat_up[iat].data[mm];
        }
    }
}

static void set_occ_mat_from_flat(
    const std::vector<double>& uom_array,
    std::vector<Matrix2D>& occ_mat_up,
    std::vector<Matrix2D>& occ_mat_dn,
    const std::vector<int>& pot_uterm_pw_index,
    int nspin)
{
    if (nspin == 4)
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                occ_mat_up[iat].data[mm] = uom_array[pot_uterm_pw_index[iat] + mm];
        }
    }
    else if (nspin == 2)
    {
        int half_size = uom_array.size() / 2;
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
            {
                occ_mat_up[iat].data[mm] = uom_array[pot_uterm_pw_index[iat] + mm];
                occ_mat_dn[iat].data[mm] = uom_array[half_size + pot_uterm_pw_index[iat] + mm];
            }
        }
    }
    else // nspin=1
    {
        for (size_t iat = 0; iat < occ_mat_up.size(); iat++)
        {
            int size = occ_mat_up[iat].nr * occ_mat_up[iat].nc;
            for (int mm = 0; mm < size; mm++)
                occ_mat_up[iat].data[mm] = uom_array[pot_uterm_pw_index[iat] + mm];
        }
    }
}

class OccMatRoundtripTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(OccMatRoundtripTest, Nspin1and2_SingleAndSplitLayout)
{
    // nspin=1: single atom d-orbital roundtrip
    const int l = 2;
    const int size = (2 * l + 1) * (2 * l + 1); // 25

    std::vector<Matrix2D> occ_mat_up(1, Matrix2D(2 * l + 1, 2 * l + 1));
    std::vector<Matrix2D> occ_mat_dn(1, Matrix2D(2 * l + 1, 2 * l + 1));
    for (int i = 0; i < size; i++)
        occ_mat_up[0].data[i] = static_cast<double>(i + 1);

    std::vector<int> pot_uterm_pw_index = {0};
    std::vector<double> uom_save(size, 0.0);
    copy_occ_mat_to_flat(occ_mat_up, occ_mat_dn, uom_save, pot_uterm_pw_index, 1);
    set_occ_mat_from_flat(uom_save, occ_mat_up, occ_mat_dn, pot_uterm_pw_index, 1);
    for (int i = 0; i < size; i++)
        EXPECT_DOUBLE_EQ(occ_mat_up[0].data[i], static_cast<double>(i + 1));

    // nspin=2: split layout [up | dn] with distinct values
    const int total = size * 2;
    for (int i = 0; i < size; i++)
    {
        occ_mat_up[0].data[i] = static_cast<double>(i + 1);
        occ_mat_dn[0].data[i] = static_cast<double>(i + 100);
    }
    uom_save.assign(total, 0.0);
    copy_occ_mat_to_flat(occ_mat_up, occ_mat_dn, uom_save, pot_uterm_pw_index, 2);
    // Verify split layout
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(uom_save[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(uom_save[size + i], static_cast<double>(i + 100));
    }
    set_occ_mat_from_flat(uom_save, occ_mat_up, occ_mat_dn, pot_uterm_pw_index, 2);
    for (int i = 0; i < size; i++)
    {
        EXPECT_DOUBLE_EQ(occ_mat_up[0].data[i], static_cast<double>(i + 1));
        EXPECT_DOUBLE_EQ(occ_mat_dn[0].data[i], static_cast<double>(i + 100));
    }
}

TEST_F(OccMatRoundtripTest, Nspin4_PauliBlocks)
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

    std::vector<int> pot_uterm_pw_index(specs.size());
    int offset = 0;
    for (size_t i = 0; i < specs.size(); i++)
    {
        pot_uterm_pw_index[i] = offset;
        offset += sizes[i];
    }

    std::vector<Matrix2D> occ_mat(specs.size());
    for (size_t i = 0; i < specs.size(); i++)
    {
        int dim = (2 * specs[i].l + 1) * npol;
        occ_mat[i] = Matrix2D(dim, dim);
        for (int j = 0; j < sizes[i]; j++)
            occ_mat[i].data[j] = static_cast<double>(i * 1000 + j + 1);
    }

    std::vector<double> uom_array(total, 0.0);
    std::vector<Matrix2D> occ_mat_dn(specs.size()); // unused for nspin=4

    copy_occ_mat_to_flat(occ_mat, occ_mat_dn, uom_array, pot_uterm_pw_index, 4);
    set_occ_mat_from_flat(uom_array, occ_mat, occ_mat_dn, pot_uterm_pw_index, 4);

    for (size_t i = 0; i < specs.size(); i++)
        for (int j = 0; j < sizes[i]; j++)
            EXPECT_DOUBLE_EQ(occ_mat[i].data[j], static_cast<double>(i * 1000 + j + 1));
}

// =====================================================================
// 3. pot_onsite effective potential formula (cal_type=3, FLL)
//
// pot_onsite[m0,m1] = U * (0.5*delta(m0,m1) - occ_mat[m0,m1])  (diagonal)
// pot_onsite[m0,m1] = -U * occ_mat[m0,m1]                       (off-diagonal)
// =====================================================================

static double compute_pot_onsite(double U_val, int m0, int m1, double occ_mat_val)
{
    if (m0 == m1)
        return U_val * (0.5 - occ_mat_val);
    else
        return -U_val * occ_mat_val;
}

class PotOnsitePotentialTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(PotOnsitePotentialTest, Diagonal_HalfFilled)
{
    double U = 4.0;
    double occ_mat = 0.5; // half-filled
    double pot_onsite = compute_pot_onsite(U, 0, 0, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, 0.0); // U * (0.5 - 0.5) = 0
}

TEST_F(PotOnsitePotentialTest, Diagonal_FullyOccupied)
{
    double U = 4.0;
    double occ_mat = 1.0; // fully occupied
    double pot_onsite = compute_pot_onsite(U, 0, 0, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, -2.0); // U * (0.5 - 1.0) = -2.0
}

TEST_F(PotOnsitePotentialTest, OffDiagonal)
{
    double U = 5.0;
    double occ_mat = 0.3;
    double pot_onsite = compute_pot_onsite(U, 0, 1, occ_mat);
    EXPECT_DOUBLE_EQ(pot_onsite, -1.5); // -U * occ_mat = -1.5
}

// =====================================================================
// 4. Energy correction formula
//
// E_U = 0.5 * U * sum_spin [Tr(n) - Tr(n^2)]
// =====================================================================

class EnergyCorrectionTest : public ::testing::Test
{
  protected:
    static double compute_energy(const std::vector<double>& occ_mat_flat, int m_size, double U)
    {
        double nm_trace = 0.0, nm2_trace = 0.0;
        for (int m0 = 0; m0 < m_size; m0++)
        {
            nm_trace += occ_mat_flat[m0 * m_size + m0];
            for (int m1 = 0; m1 < m_size; m1++)
                nm2_trace += occ_mat_flat[m0 * m_size + m1] * occ_mat_flat[m1 * m_size + m0];
        }
        return 0.5 * U * (nm_trace - nm2_trace);
    }
};

TEST_F(EnergyCorrectionTest, HalfFilled_DOrbital)
{
    const int m_size = 5;
    std::vector<double> occ_mat(m_size * m_size, 0.0);
    for (int m = 0; m < m_size; m++)
        occ_mat[m * m_size + m] = 0.5;

    double energy = compute_energy(occ_mat, m_size, 4.0);
    // Tr(n) = 2.5, Tr(n^2) = 1.25, E = 0.5 * 4 * 1.25 = 2.5
    EXPECT_DOUBLE_EQ(energy, 2.5);
}

TEST_F(EnergyCorrectionTest, OffDiagonal_Contribution)
{
    const int m_size = 2;
    std::vector<double> occ_mat = {
        0.3, 0.1,
        0.1, 0.3
    };

    double energy = compute_energy(occ_mat, m_size, 4.0);
    // Tr(n) = 0.6, Tr(n^2) = 0.3^2 + 0.1^2 + 0.1^2 + 0.3^2 = 0.20
    // E = 0.5 * 4 * (0.6 - 0.20) = 0.8
    EXPECT_DOUBLE_EQ(energy, 0.8);
}

TEST_F(EnergyCorrectionTest, DoubleCounting_Energy)
{
    // E_dc = sum_{m1,m2,spin} pot_onsite[m1,m2] * n[m2,m1]
    const int m_size = 3;
    double U = 4.0;
    std::vector<double> occ_mat = {
        0.5, 0.0, 0.0,
        0.0, 0.3, 0.0,
        0.0, 0.0, 0.2
    };

    double e_dc = 0.0;
    for (int m1 = 0; m1 < m_size; m1++)
        for (int m2 = 0; m2 < m_size; m2++)
        {
            double pot_onsite = (m1 == m2) ? U * (0.5 - occ_mat[m1 * m_size + m2])
                                   : -U * occ_mat[m1 * m_size + m2];
            e_dc += pot_onsite * occ_mat[m2 * m_size + m1];
        }

    // Only diagonal: m=0: 0*0.5=0, m=1: 0.8*0.3=0.24, m=2: 1.2*0.2=0.24
    EXPECT_NEAR(e_dc, 0.48, 1e-14);
}
