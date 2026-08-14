#include <cmath>
#include <complex>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_basis/module_ao/parallel_orbitals.h"

/************************************************
 * T3: closed Frobenius trace equivalence
 *
 * For a Hermitian DMK(k), a real DMR obtained via cal_DMR (e^{-ikR}), and a
 * paired Hermitian operator O(R) (O(iat2,iat1,-R) = O(iat1,iat2,R)^\dagger),
 * the real-space contraction
 *
 *     W_R = \sum_{mu,nu,R} D(mu,nu;R) * O(mu,nu;R)
 *
 * must equal the k-space contraction
 *
 *     W_k = \sum_k Re Tr( DMK(k) * O(k) ),   O(k) = folding_HR(O(R), k).
 *
 * This holds even when O(R) is NOT symmetric within the same R (the key case
 * for overlap/kinetic derivative operators): the protection comes from the
 * closed trace + per-k Hermiticity of O(k), not from same-R symmetry.
 * See docs/dm_dmk_dmr_action_plan.md (T3) and
 * docs/nao_lcao_force_stress_derivation.md Sec. 5.
 ************************************************/

namespace
{
// value builders for T = double or std::complex<double>
template <typename T>
T mk_value(const double re, const double im);
template <>
double mk_value<double>(const double re, const double)
{
    return re;
}
template <>
std::complex<double> mk_value<std::complex<double>>(const double re, const double im)
{
    return std::complex<double>(re, im);
}

// conjugate transpose helper: identity for double, conj for complex
template <typename T>
T conj_val(const T& v);
template <>
double conj_val<double>(const double& v)
{
    return v;
}
template <>
std::complex<double> conj_val<std::complex<double>>(const std::complex<double>& v)
{
    return std::conj(v);
}

// Hermitian matrix H (row-major, n*n), deterministic pseudo-random values
void fill_hermitian(std::vector<std::complex<double>>& H, const int n, const unsigned seed)
{
    H.assign(n * n, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            const double re = 0.5 * std::sin(double(i * 7 + j * 13 + seed * 3));
            const double im = (i == j) ? 0.0 : 0.5 * std::cos(double(i * 11 + j * 5 + seed * 7));
            H[i * n + j] = std::complex<double>(re, im);
            H[j * n + i] = std::complex<double>(re, -im);
        }
    }
}

// HContainer with (iat1,iat2,R); if paired, also (iat2,iat1,-R)
template <typename T>
hamilt::HContainer<T> build_r_container(
    const Parallel_Orbitals* paraV,
    const std::vector<std::pair<int, int>>& pairs,
    const std::vector<ModuleBase::Vector3<int>>& r_list,
    const bool paired)
{
    hamilt::HContainer<T> hR(paraV);
    for (const auto& pr: pairs)
    {
        for (const auto& R: r_list)
        {
            hR.insert_pair(hamilt::AtomPair<T>(pr.first, pr.second, R, paraV));
            if (paired)
            {
                hR.insert_pair(hamilt::AtomPair<T>(pr.second, pr.first, -R, paraV));
            }
        }
    }
    hR.allocate(nullptr, true);
    return hR;
}

// Build a paired Hermitian operator O(R) on the container:
// - forward blocks (iat1,iat2,R) filled with deterministic values,
//   symmetric within R if symmetric_in_R is true,
// - reverse blocks (iat2,iat1,-R) set to the transpose (real) or conjugate
//   transpose (complex) of the forward blocks.
template <typename T>
hamilt::HContainer<T> build_operator_container(
    const Parallel_Orbitals* paraV,
    const std::vector<std::pair<int, int>>& pairs,
    const std::vector<ModuleBase::Vector3<int>>& r_list,
    const unsigned seed,
    const bool symmetric_in_R)
{
    auto hR = build_r_container<T>(paraV, pairs, r_list, false);
    for (int iap = 0; iap < hR.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<T>& ap = hR.get_atom_pair(iap);
        const int row_size = ap.get_row_size();
        const int col_size = ap.get_col_size();
        const bool self_pair = (ap.get_atom_i() == ap.get_atom_j());
        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
            T* mat = ap.get_pointer(ir);
            for (int i = 0; i < row_size; ++i)
            {
                for (int j = 0; j < col_size; ++j)
                {
                    double re = 0.3 * std::sin(double(i * 3 + j * 5 + iap * 7 + ir * 11 + seed));
                    double im = 0.3 * std::cos(double(i * 5 + j * 3 + iap * 11 + ir * 7 + seed));
                    mat[i * col_size + j] = mk_value<T>(re, im);
                    // the R=(0,0,0) self block must be symmetric (real) /
                    // Hermitian (complex): its reverse partner is itself
                    const bool force_hermitian = (self_pair && R.x == 0 && R.y == 0 && R.z == 0);
                    if (symmetric_in_R || force_hermitian)
                    {
                        mat[j * col_size + i] = conj_val<T>(mat[i * col_size + j]);
                    }
                }
            }
        }
    }
    // insert reverse blocks O(iat2,iat1,-R) = O(iat1,iat2,R)^T (real) or ^\dagger
    for (const auto& pr: pairs)
    {
        for (const auto& R: r_list)
        {
            // the R=(0,0,0) self block is its own reverse partner
            if (pr.first == pr.second && R.x == 0 && R.y == 0 && R.z == 0)
            {
                continue;
            }
            const hamilt::BaseMatrix<T>* fwd = hR.find_matrix(pr.first, pr.second, R);
            EXPECT_NE(fwd, nullptr) << "forward block missing for (" << pr.first << "," << pr.second << "," << R.x << "," << R.y << "," << R.z << ")";
            if (fwd == nullptr)
            {
                continue;
            }
            const int row_size = fwd->get_row_size();
            const int col_size = fwd->get_col_size();
            hamilt::AtomPair<T> ap_rev(pr.second, pr.first, -R, paraV);
            hamilt::BaseMatrix<T>& rev_mat = ap_rev.get_HR_values(-R.x, -R.y, -R.z);
            rev_mat.allocate(nullptr, true);
            T* rev = rev_mat.get_pointer();
            for (int i = 0; i < row_size; ++i)
            {
                for (int j = 0; j < col_size; ++j)
                {
                    const T fwd_val = fwd->get_pointer()[i * col_size + j];
                    rev[j * col_size + i] = conj_val<T>(fwd_val);
                }
            }
            hR.insert_pair(ap_rev);
        }
    }
    hR.allocate(nullptr, false);
    return hR;
}

// W_R = sum over container (pairs, R) of sum_{mu,nu} D(mu,nu;R) * O(mu,nu;R)
template <typename T>
std::complex<double> compute_WR(const hamilt::HContainer<double>& dmr, const hamilt::HContainer<T>& hr)
{
    std::complex<double> wr(0.0, 0.0);
    for (int iap = 0; iap < dmr.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& ap = dmr.get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        const int row_size = ap.get_row_size();
        const int col_size = ap.get_col_size();
        for (int ir = 0; ir < ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R = ap.get_R_index(ir);
            const hamilt::BaseMatrix<double>* dmat = ap.find_matrix(R);
            const hamilt::BaseMatrix<T>* omat = hr.find_matrix(iat1, iat2, R);
            EXPECT_NE(omat, nullptr) << "O missing for (" << iat1 << "," << iat2 << "," << R.x << "," << R.y << "," << R.z << ")";
            if (omat == nullptr)
            {
                continue;
            }
            for (int i = 0; i < row_size; ++i)
            {
                for (int j = 0; j < col_size; ++j)
                {
                    wr += dmat->get_pointer()[i * col_size + j] * omat->get_pointer()[i * col_size + j];
                }
            }
        }
    }
    return wr;
}

// W_k = sum_k Re Tr( DMK(k) * folding_HR(O(R), k) )
template <typename T>
double compute_Wk(const std::vector<std::vector<std::complex<double>>>& M,
                  const int nk,
                  const hamilt::HContainer<T>& hr,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                  const int ntot)
{
    double wk = 0.0;
    for (int ik = 0; ik < nk; ++ik)
    {
        std::vector<std::complex<double>> ok(ntot * ntot, std::complex<double>(0.0, 0.0));
        hamilt::folding_HR(hr, ok.data(), kvec_d[ik], ntot, 0);
        std::complex<double> tr(0.0, 0.0);
        for (int mu = 0; mu < ntot; ++mu)
        {
            for (int nu = 0; nu < ntot; ++nu)
            {
                tr += M[ik][mu * ntot + nu] * ok[nu * ntot + mu];
            }
        }
        // Tr(DMK(k) O(k)) is real for Hermitian DMK(k) and Hermitian O(k)
        EXPECT_NEAR(tr.imag(), 0.0, 1e-10) << "trace not real at ik=" << ik;
        wk += tr.real();
    }
    return wk;
}
} // namespace

class DMTraceTest : public ::testing::Test
{
  protected:
    enum : int
    {
        nat = 2,
        nw = 3, // orbitals per atom
        ntot = nat * nw
    };

    Parallel_Orbitals* paraV;

    void SetUp() override
    {
        paraV = new Parallel_Orbitals();
        const int iat2iwt[nat] = {0, nw};
        paraV->init(ntot, ntot, 2, MPI_COMM_WORLD);
        paraV->set_atomic_trace(iat2iwt, nat, ntot);
    }

    void TearDown() override
    {
        delete paraV;
    }
};

TEST_F(DMTraceTest, T3_closed_trace_equivalence)
{
    const int nk = 3;
    const std::vector<ModuleBase::Vector3<double>> kvec_d = {
        ModuleBase::Vector3<double>(0.1, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.2, 0.2, 0.0),
        ModuleBase::Vector3<double>(0.3, 0.1, 0.2)};

    // Hermitian DMK per k
    std::vector<std::vector<std::complex<double>>> M(nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        fill_hermitian(M[ik], ntot, ik);
    }

    elecstate::DensityMatrix<std::complex<double>, double> DM(paraV, 1, kvec_d, nk);
    for (int ik = 0; ik < nk; ++ik)
    {
        for (int i = 0; i < ntot; ++i)
        {
            for (int j = 0; j < ntot; ++j)
            {
                // set_DMK stores the transposed element (see density_matrix.h)
                DM.set_DMK(1, ik, j, i, M[ik][i * ntot + j]);
            }
        }
    }

    // paired DMR container: cross pair (0,1) and self pairs
    const std::vector<std::pair<int, int>> pairs = {{0, 1}, {0, 0}, {1, 1}};
    const std::vector<ModuleBase::Vector3<int>> r_list = {
        ModuleBase::Vector3<int>(0, 0, 0),
        ModuleBase::Vector3<int>(1, 0, 0),
        ModuleBase::Vector3<int>(0, 1, 0),
        ModuleBase::Vector3<int>(0, 0, 1)};

    hamilt::HContainer<double> dmr = build_r_container<double>(paraV, pairs, r_list, true);
    DM.init_DMR(dmr);
    DM.cal_DMR();

    // case 1: real symmetric O(R) (overlap-like)
    {
        auto hr = build_operator_container<double>(paraV, pairs, r_list, 1, true);
        const std::complex<double> wr = compute_WR(dmr, hr);
        EXPECT_NEAR(wr.imag(), 0.0, 1e-10) << "case 1: W_R should be real";
        const double wk = compute_Wk(M, nk, hr, kvec_d, ntot);
        EXPECT_NEAR(wr.real(), wk, 1e-10) << "case 1: symmetric real O(R)";
    }
    // case 2: non-symmetric real O(R) within the same R (derivative-like).
    // The equivalence must still hold: protection comes from the closed trace
    // + per-k Hermiticity of O(k), NOT from same-R symmetry.
    {
        auto hr = build_operator_container<double>(paraV, pairs, r_list, 2, false);
        const std::complex<double> wr = compute_WR(dmr, hr);
        EXPECT_NEAR(wr.imag(), 0.0, 1e-10) << "case 2: W_R should be real";
        const double wk = compute_Wk(M, nk, hr, kvec_d, ntot);
        EXPECT_NEAR(wr.real(), wk, 1e-10) << "case 2: non-symmetric real O(R)";
    }
    // case 3: complex Hermitian O(R) (multi-k)
    {
        auto hr = build_operator_container<std::complex<double>>(paraV, pairs, r_list, 3, false);
        const std::complex<double> wr = compute_WR(dmr, hr);
        EXPECT_NEAR(wr.imag(), 0.0, 1e-10) << "case 3: W_R should be real";
        const double wk = compute_Wk(M, nk, hr, kvec_d, ntot);
        EXPECT_NEAR(wr.real(), wk, 1e-10) << "case 3: complex Hermitian O(R)";
    }
}

// Stub required by the DensityMatrix member functions instantiated from
// density_matrix_io.cpp (Record_adj::cal_adj is not linked into this
// unit-test target; see also tmp_mocks.cpp).
#include "source_lcao/record_adj.h"
Record_adj::Record_adj() {}
Record_adj::~Record_adj() {}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
