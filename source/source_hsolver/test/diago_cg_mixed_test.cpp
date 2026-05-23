#include "gtest/gtest.h"
#include "source_hsolver/diago_cg.h"
#include <complex>
#include <random>
#include <vector>

using Complex = std::complex<double>;

static void make_hermitian(int n, std::vector<Complex>& H)
{
    H.resize(n * n);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double real = dist(rng);
            const double imag = (i == j ? 0.0 : dist(rng));
            H[i * n + j] = Complex(real, imag);
            H[j * n + i] = std::conj(H[i * n + j]);
        }
    }
}

static void make_random_psi(int nband, int dim, std::vector<Complex>& psi)
{
    psi.resize(static_cast<size_t>(nband) * dim);
    std::mt19937_64 rng(54321);
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    for (int i = 0; i < nband * dim; ++i) {
        psi[i] = Complex(dist(rng), dist(rng));
    }
}

static void apply_hamiltonian(const std::vector<Complex>& H,
                              int n,
                              const Complex* psi_in,
                              Complex* hpsi_out,
                              int ld, int nvec)
{
    for (int v = 0; v < nvec; ++v) {
        const Complex* psi_vec = psi_in + static_cast<size_t>(v) * ld;
        Complex* out_vec = hpsi_out + static_cast<size_t>(v) * ld;
        for (int i = 0; i < n; ++i) {
            Complex sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += H[static_cast<size_t>(i) * n + j] * psi_vec[j];
            }
            out_vec[i] = sum;
        }
    }
}

static void apply_overlap(const Complex* psi_in,
                          Complex* spsi_out,
                          int ld,
                          int nvec)
{
    for (int i = 0; i < nvec * ld; ++i) {
        spsi_out[i] = psi_in[i];
    }
}

TEST(DiagoCGMixedTest, MixedPrecisionMatchesDouble)
{
    const int dim = 8;
    const int nband = 3;
    const int ld_psi = dim;

    std::vector<Complex> H;
    make_hermitian(dim, H);

    std::vector<Complex> psi_initial;
    make_random_psi(nband, dim, psi_initial);

    std::vector<Complex> psi_double = psi_initial;
    std::vector<Complex> psi_mixed = psi_initial;
    std::vector<double> eigen_double(nband, 0.0);
    std::vector<double> eigen_mixed(nband, 0.0);

    auto hpsi_func = [&H, dim](Complex* psi_in, Complex* hpsi_out, const int ld, const int nvec) {
        apply_hamiltonian(H, dim, psi_in, hpsi_out, ld, nvec);
    };
    auto spsi_func = [](Complex* psi_in, Complex* spsi_out, const int ld, const int nvec) {
        apply_overlap(psi_in, spsi_out, ld, nvec);
    };

    std::vector<double> ethr_band(nband, 1e-6);

    hsolver::DiagoCG<Complex> cg_double(
        "pw",
        "nscf",
        false,
        hsolver::DiagoCG<Complex>::SubspaceFunc(),
        1e-6,
        200,
        1);
    cg_double.diag(hpsi_func,
                   spsi_func,
                   ld_psi,
                   nband,
                   dim,
                   psi_double.data(),
                   eigen_double.data(),
                   ethr_band,
                   nullptr);

    hsolver::DiagoCG<Complex> cg_mixed(
        "pw",
        "nscf",
        false,
        hsolver::DiagoCG<Complex>::SubspaceFunc(),
        1e-6,
        200,
        1,
        hsolver::PrecisionMode::kMixed);
    cg_mixed.diag(hpsi_func,
                  spsi_func,
                  ld_psi,
                  nband,
                  dim,
                  psi_mixed.data(),
                  eigen_mixed.data(),
                  ethr_band,
                  nullptr);

    for (int i = 0; i < nband; ++i) {
        EXPECT_NEAR(eigen_double[i], eigen_mixed[i], 1e-6) << "Index=" << i;
    }
}
