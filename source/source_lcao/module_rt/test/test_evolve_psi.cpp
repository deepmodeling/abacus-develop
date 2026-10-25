#include "source_lcao/module_rt/evolve_psi.h"
#include "source_base/module_external/scalapack_connector.h"
#include "tddft_test.h"

#include <algorithm>
#include <array>
#include <complex>
#include <fstream>

namespace
{
using Complex = std::complex<double>;

// A two-orbital Hamiltonian with changing overlap and zero energy. This keeps
// the normalization oracle independent of the propagator's phase convention.
class EndpointHamiltonian : public hamilt::Hamilt<Complex>
{
  public:
    std::array<Complex, 4> h = {{0.0, 0.0, 0.0, 0.0}};
    std::array<Complex, 4> s = {{2.0, Complex(0.0, -0.4), Complex(0.0, 0.4), 1.0}};

    void matrix(hamilt::MatrixBlock<Complex>& hk, hamilt::MatrixBlock<Complex>& sk) override
    {
        hk.p = h.data();
        sk.p = s.data();
        hk.row = sk.row = 2;
        hk.col = sk.col = 2;
    }
};

class EvolvePsiTest : public testing::Test
{
  protected:
    Parallel_Orbitals pv;
    EndpointHamiltonian hamiltonian;
    std::array<Complex, 4> previous_s = {{1.0, 0.0, 0.0, 1.0}};
    std::array<Complex, 2> previous_psi = {{1.0, Complex(0.0, 0.5)}};
    const module_rt::EijOutputInfo output_info{0.0, 0, 0, 0, 0, 0};
    std::ofstream log;

    void SetUp() override
    {
        const int nlocal = 2;
        const int nband = 1;
        const int block = 1;
        const int zero = 0;
        const int lld = 2;
        int info = 0;
        pv.nloc = 4;
        pv.nloc_wfc = 2;
        pv.nrow = pv.ncol = 2;
        pv.nrow_bands = pv.ncol_bands = 1;
        pv.dim0 = pv.dim1 = pv.nb = 1;
        pv.blacs_ctxt = ictxt;
        pv.set_coord(0, 0);
        descinit_(pv.desc, &nlocal, &nlocal, &block, &block, &zero, &zero, &ictxt, &lld, &info);
        ASSERT_EQ(info, 0);
        descinit_(pv.desc_wfc, &nlocal, &nband, &block, &block, &zero, &zero, &ictxt, &lld, &info);
        ASSERT_EQ(info, 0);
        descinit_(pv.desc_Eij, &nband, &nband, &block, &block, &zero, &zero, &ictxt, &lld, &info);
        ASSERT_EQ(info, 0);
    }

    void expect_endpoint_norm(const Complex* psi)
    {
        const auto& s = hamiltonian.s;
        const Complex norm = std::conj(psi[0]) * (s[0] * psi[0] + s[2] * psi[1])
                           + std::conj(psi[1]) * (s[1] * psi[0] + s[3] * psi[1]);
        EXPECT_NEAR(norm.real(), 1.0, 1e-12);
        EXPECT_NEAR(norm.imag(), 0.0, 1e-12);
    }

    void run_tensor_lapack()
    {
        ct::Tensor psi(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({1, 2}));
        ct::Tensor last_psi(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({1, 2}));
        ct::Tensor last_h(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({4}));
        ct::Tensor last_s(ct::DataType::DT_COMPLEX_DOUBLE, ct::TensorShape({4}));
        ct::Tensor energy(ct::DataType::DT_DOUBLE, ct::TensorShape({1}));
        last_h.zero();
        std::copy(previous_psi.begin(), previous_psi.end(), last_psi.data<Complex>());
        std::copy(previous_s.begin(), previous_s.end(), last_s.data<Complex>());
        CublasMpResources resources;
        const bool use_lapack = true;
        module_rt::evolve_psi_tensor<base_device::DEVICE_CPU>(
            1, 2, &pv, &hamiltonian, psi, last_psi, last_h, last_s, energy, 0, log, 0,
            use_lapack, resources, output_info);
        expect_endpoint_norm(psi.data<Complex>());
        EXPECT_NEAR(energy.data<double>()[0], 0.0, 1e-12);
    }
};

TEST_F(EvolvePsiTest, CpuUsesEndpointOverlap)
{
    std::array<Complex, 2> psi;
    std::array<Complex, 4> last_h = {{0.0, 0.0, 0.0, 0.0}};
    double energy = 0.0;
    module_rt::evolve_psi(1, 2, &pv, &hamiltonian, psi.data(), previous_psi.data(),
                          last_h.data(), previous_s.data(), nullptr, false, &energy, 0, log, 0, output_info);
    expect_endpoint_norm(psi.data());
    EXPECT_NEAR(energy, 0.0, 1e-12);
}

TEST_F(EvolvePsiTest, LapackUsesEndpointOverlap)
{
    run_tensor_lapack();
}
} // namespace
