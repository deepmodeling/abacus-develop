#include "../potentials/xc_kernel.h"
#include "../utils/lr_util.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_output/cube_io.h"
#include <gtest/gtest.h>
#include "source_base/parallel_cell.h"
#include <limits>

// Keep the actual LibXC evaluation and coefficient construction. Geometry ownership,
// FFTs and grid masks are isolated; no local grid point may be accessed on an empty rank.
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}
Charge::Charge() {}
Charge::~Charge() {}

namespace ModuleBase
{
FFT_Bundle::~FFT_Bundle() {}
} // namespace ModuleBase

namespace ModulePW
{
PW_Basis::PW_Basis() {}
PW_Basis::~PW_Basis() {}
void PW_Basis::initgrids(double, ModuleBase::Matrix3, double) {}
void PW_Basis::initgrids(
    double,
    ModuleBase::Matrix3,
    int,
    int,
    int)
{
}
void PW_Basis::distribute_r() {}
template <typename T>
void PW_Basis::real2recip(
    const T *in,
    std::complex<T> *out,
    const bool,
    const T) const
{
    for (int point = 0; point < npw; ++point)
    {
        out[point] = in[point];
    }
}
template void PW_Basis::real2recip<double>(
    const double *,
    std::complex<double> *,
    bool,
    double) const;
template <typename T>
void PW_Basis::recip2real(
    const std::complex<T> *in,
    std::complex<T> *out,
    const bool,
    const T) const
{
    for (int point = 0; point < nrxx; ++point)
    {
        out[point] = in[point];
    }
}
template void PW_Basis::recip2real<double>(
    const std::complex<double> *,
    std::complex<double> *,
    bool,
    double) const;
} // namespace ModulePW

void XC_Functional::grad_rho(
    const std::complex<double> *,
    ModuleBase::Vector3<double> *out,
    const ModulePW::PW_Basis *basis,
    const double tpiba)
{
    for (int point = 0; point < basis->nrxx; ++point)
    {
        out[point] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3) * tpiba;
    }
}
void XC_Functional::grad_dot(
    const ModuleBase::Vector3<double> *,
    double *out,
    const ModulePW::PW_Basis *basis,
    const double)
{
    for (int point = 0; point < basis->nrxx; ++point)
    {
        out[point] = 0.0;
    }
}
namespace ModuleIO
{
bool read_vdata_palgrid(
    const Parallel_Grid &,
    const int,
    std::ostream &,
    const std::string &,
    double *const,
    const int)
{
    throw std::logic_error("This coefficient test does not read kernel files");
}
} // namespace ModuleIO
namespace XC_Functional_Libxc
{
std::vector<double> cal_sgn(
    const double,
    const double,
    const xc_func_type &,
    const int nspin,
    const std::size_t nrxx,
    const std::vector<double> &,
    const std::vector<double> &)
{
    return std::vector<double>(nrxx * nspin, 1.0);
}
} // namespace XC_Functional_Libxc

/** @brief Exercise coefficient availability with real LibXC on zero or one local grid point. */
class KernelXCTest : public testing::Test
{
  protected:
    ModulePW::PW_Basis basis;
    UnitCell cell;
    Charge charge;
    Parallel_Grid grid;
    double density_value = 0.3;
    double *density[2] = {&density_value, &density_value};

    void SetUp() override
    {
        basis.nrxx = 0;
        basis.npw = 0;
        basis.nmaxgr = 0;
        cell.omega = 1.0;
        cell.tpiba = 1.0;
        charge.nrxx = 0;
        charge.rho = density;
    }

    std::shared_ptr<const LR::KernelXC> build(const std::string &name, const int spins, const int requested)
    {
        const std::vector<std::string> initialization = {"default"};
        return LR::build_xc_kernel(
            name,
            basis,
            cell,
            charge,
            grid,
            spins,
            false,
            requested,
            initialization);
    }

    /**
     * @brief Respect the installed LibXC HSE derivative capability in coefficient tests.
     * @param[in] spins Native spin-component count, one or two.
     * @return Whether every HSE component advertises third derivatives.
     */
    bool hse_has_kxc(const int spins)
    {
        XC_Functional::set_xc_type("hse");
        const int polarization = spins == 1 ? XC_UNPOLARIZED : XC_POLARIZED;
        const double alpha = XC_Functional::get_hybrid_alpha();
        const double omega = XC_Functional::get_hse_omega();
        auto functions = XC_Functional_Libxc::init_func(
            XC_Functional::get_func_id(),
            polarization,
            alpha,
            omega);
        bool available = true;
        for (const auto &functional : functions)
        {
            available = available && (functional.info->flags & XC_FLAGS_HAVE_KXC) != 0;
        }
        XC_Functional_Libxc::finish_func(functions);
        return available;
    }
};

TEST_F(KernelXCTest, HseGroundResponseUsesUnmodifiedSigma)
{
    basis.nrxx = 1;
    basis.npw = 1;
    basis.nmaxgr = 1;
    charge.nrxx = 1;
    cell.tpiba = 0.001;
    const std::vector<std::string> initialization = {"default"};
    const auto ground = LR::build_gs_xc_kernel(
        "hse",
        basis,
        cell,
        charge,
        grid,
        1,
        false,
        initialization);
    const auto response = build("hse", 1, LR::KernelXC::NoGxc);
    const auto gradient = ModuleBase::Vector3<double>(0.1, 0.2, 0.3) * cell.tpiba;
    const double sigma = gradient * gradient;
    const double alpha = XC_Functional::get_hybrid_alpha();
    const double omega = XC_Functional::get_hse_omega();
    auto functions = XC_Functional_Libxc::init_func(
        XC_Functional::get_func_id(),
        XC_UNPOLARIZED,
        alpha,
        omega);
    ASSERT_EQ(functions.size(), 1);
    double frr = 0.0;
    double frs = 0.0;
    double fss = 0.0;
    xc_func_set_dens_threshold(&functions[0], 1e-14);
    xc_gga_fxc(
        &functions[0],
        1,
        &density_value,
        &sigma,
        &frr,
        &frs,
        &fss);
    XC_Functional_Libxc::finish_func(functions);
    EXPECT_FALSE(ground->gs_sigma_fallback_used());
    EXPECT_DOUBLE_EQ(ground->v2rho2[0], frr);
    EXPECT_DOUBLE_EQ(ground->v2rhosigma[0], frs);
    EXPECT_DOUBLE_EQ(ground->v2sigma2[0], fss);
    EXPECT_NE(ground->v2rhosigma[0], response->v2rhosigma[0]);
    EXPECT_THROW(ground->gxc(false), std::runtime_error);
    EXPECT_THROW(
        LR::KernelXC(
            LR::KernelXC::Purpose::GroundStateResponse,
            basis,
            cell,
            charge,
            grid,
            1,
            "hse",
            initialization,
            false,
            LR::KernelXC::Singlet),
        std::invalid_argument);
}

// A genuine LibXC HSE small-gradient failure. Only one rank has points; all ranks
// must retry together, including empty ranks. Compare every spin coefficient with
// a separately evaluated floored LR kernel, including finite raw entries.
TEST_F(KernelXCTest, HseGroundResponseFallbackIncludesEmptyRanksAndAllSpinBlocks)
{
    const auto domain = ModuleBase::world_comm_domain();
    const int candidate_rank = domain.rank();
    for (const int spins : {1, 2})
    {
        const bool owns_points = candidate_rank == 0;
        basis.nrxx = owns_points ? 1 : 0;
        basis.npw = basis.nrxx;
        basis.nmaxgr = basis.nrxx;
        charge.nrxx = basis.nrxx;
        cell.tpiba = 0.; // The fixture's mock gradient is zero: sigma = 0.
        const std::vector<std::string> initialization = {"default"};
        const auto ground = LR::build_gs_xc_kernel(
            "hse", basis, cell, charge, grid, spins, false, initialization);
        const auto floored = build("hse", spins, LR::KernelXC::NoGxc);
        // A library with a corrected sigma=0 evaluation need not take the fallback.
        // The installed KXC-enabled LibXC 5.2.3 exercises it for both spin layouts.
        if (!ground->gs_sigma_fallback_used())
        {
            GTEST_SKIP() << "Installed LibXC evaluates this raw HSE point finitely";
        }
        EXPECT_EQ(ground->vsigma, floored->vsigma);
        EXPECT_EQ(ground->v2rho2, floored->v2rho2);
        EXPECT_EQ(ground->v2rhosigma, floored->v2rhosigma);
        EXPECT_EQ(ground->v2sigma2, floored->v2sigma2);
        EXPECT_TRUE(ground->sigma_jacobian().empty());
        EXPECT_THROW(ground->gxc(false), std::runtime_error);
    }
}

TEST_F(KernelXCTest, HseGroundResponseRetryFailsCollectivelyIfStillNonfinite)
{
    const auto domain = ModuleBase::world_comm_domain();
    const bool owns_points = domain.rank() == 0;
    basis.nrxx = owns_points ? 1 : 0;
    basis.npw = basis.nrxx;
    basis.nmaxgr = basis.nrxx;
    charge.nrxx = basis.nrxx;
    cell.tpiba = std::numeric_limits<double>::quiet_NaN();
    const std::vector<std::string> initialization = {"default"};
    EXPECT_THROW(LR::build_gs_xc_kernel(
        "hse", basis, cell, charge, grid, 1, false, initialization), std::runtime_error);
}

TEST_F(KernelXCTest, SigmaFloorJacobianPreservesOtherFunctionals)
{
    const std::vector<double> sigma = {-2.0,
        0.5e-6,
        1e-6,
        2e-6};
    std::vector<double> buffer;
    std::vector<double> jacobian;
    const auto &original = LR_Util::prepare_xc_sigma(
        sigma,
        false,
        buffer,
        jacobian);
    EXPECT_EQ(&original, &sigma);
    EXPECT_TRUE(jacobian.empty());
    const auto &stabilized = LR_Util::prepare_xc_sigma(
        sigma,
        true,
        buffer,
        jacobian);
    EXPECT_EQ(
        stabilized, (std::vector<double>{1e-6,
                        1e-6,
                        1e-6,
                        2e-6}));
    EXPECT_EQ(
        jacobian, (std::vector<double>{0.0,
                      0.0,
                      0.0,
                      1.0}));
    const std::vector<double> empty;
    EXPECT_TRUE(
        LR_Util::prepare_xc_sigma(
            empty,
            true,
            buffer,
            jacobian)
            .empty());
    EXPECT_TRUE(jacobian.empty());
}


TEST_F(KernelXCTest, HseClippedKernelGradientMatchesIndependentFiniteDifference)
{
    if (!hse_has_kxc(1))
    {
        GTEST_SKIP() << "Installed LibXC does not advertise HSE third derivatives";
    }
    basis.nrxx = 1;
    basis.npw = 1;
    basis.nmaxgr = 1;
    charge.nrxx = 1;
    cell.tpiba = 0.001;
    const auto kernel = build("hse", 1, LR::KernelXC::Singlet);
    const auto &coefficients = kernel->gxc(false);
    const ModuleBase::Vector3<double> base_gradient(0.1, 0.2, 0.3);
    const ModuleBase::Vector3<double> perturbation_gradient(0.03, -0.02, 0.04);
    const double u = 0.4;
    const auto gradient = base_gradient * cell.tpiba;
    const double t = gradient * perturbation_gradient;
    const double q = perturbation_gradient * perturbation_gradient;
    const double radial = coefficients.e_s2[0] * u * u + coefficients.e_st[0] * u * t + coefficients.e_t2[0] * t * t +
                          coefficients.e_q[0] * q;
    const auto& second_order = kernel->closed_shell_gga(false);
    const double explicit_gradient = 8.0 * (second_order.rho_sigma[0] * u + second_order.sigma_sigma[0] * t);
    const auto derivative = gradient * radial + perturbation_gradient * explicit_gradient;
    const double analytic = (derivative * base_gradient) / 4.0;
    // Differentiate the second-order LR quadratic form, without using any third derivative.
    const auto quadratic = [&](const double scale) -> double
    {
        cell.tpiba = scale;
        const auto second = build("hse", 1, LR::KernelXC::NoGxc);
        const double projection = (base_gradient * scale) * perturbation_gradient;
        return second->v2rho2[0] * u * u + 4.0 * second->v2rhosigma[0] * u * projection +
               4.0 * second->v2sigma2[0] * projection * projection + 2.0 * second->vsigma[0] * q;
    };
    for (const double h : {1e-5, 5e-6})
    {
        const double plus = quadratic(0.001 + h);
        const double minus = quadratic(0.001 - h);
        const double finite_difference = (plus - minus) / (2.0 * h);
        EXPECT_NEAR(analytic, finite_difference, 1e-9);
    }
    EXPECT_DOUBLE_EQ(radial, 0.0);
    EXPECT_NE(explicit_gradient, 0.0);
}

TEST_F(KernelXCTest, PolarizedHseClippingMasksOnlyImplicitSigmaDerivatives)
{
    if (!hse_has_kxc(2))
    {
        GTEST_SKIP() << "Installed LibXC does not advertise polarized HSE third derivatives";
    }
    basis.nrxx = 1;
    basis.npw = 1;
    basis.nmaxgr = 1;
    charge.nrxx = 1;
    cell.tpiba = 0.001;
    const auto kernel = build("hse", 2, LR::KernelXC::BothSpins);
    for (const bool triplet : {false, true})
    {
        const auto &coefficients = kernel->gxc(triplet);
        EXPECT_DOUBLE_EQ(coefficients.e_s2[0], 0.0);
        EXPECT_DOUBLE_EQ(coefficients.e_st[0], 0.0);
        EXPECT_DOUBLE_EQ(coefficients.e_t2[0], 0.0);
        EXPECT_DOUBLE_EQ(coefficients.e_q[0], 0.0);
        EXPECT_NE(kernel->closed_shell_gga(triplet).rho_sigma[0], 0.0);
    }
    const auto pbe = build("pbe", 2, LR::KernelXC::BothSpins);
    EXPECT_NE(pbe->gxc(false).e_s2[0], 0.0);
}

TEST_F(KernelXCTest, RequestedCoefficientsExistOnEmptyRanks)
{
    for (const std::string &name : {"lda", "pbe"})
    {
        for (const int spins : {1, 2})
        {
            const auto kernel = build(name, spins, LR::KernelXC::BothSpins);
            EXPECT_TRUE(kernel->gxc(false).a_s2.empty());
            EXPECT_TRUE(kernel->gxc(true).a_s2.empty());
            EXPECT_EQ(kernel->is_gga(), name == "pbe");
        }
    }
}

TEST_F(KernelXCTest, MissingSpinStillThrowsOnEmptyRanks)
{
    const auto absent = build("pbe", 2, LR::KernelXC::NoGxc);
    EXPECT_THROW(absent->gxc(false), std::runtime_error);
    EXPECT_THROW(absent->gxc(true), std::runtime_error);
    const auto singlet = build("pbe", 2, LR::KernelXC::Singlet);
    EXPECT_NO_THROW(singlet->gxc(false));
    EXPECT_THROW(singlet->gxc(true), std::runtime_error);
    const auto triplet = build("pbe", 2, LR::KernelXC::Triplet);
    EXPECT_THROW(triplet->gxc(false), std::runtime_error);
    EXPECT_NO_THROW(triplet->gxc(true));
    const auto nonlocal = build("hf", 2, LR::KernelXC::BothSpins);
    EXPECT_THROW(nonlocal->gxc(false), std::runtime_error);
}

TEST_F(KernelXCTest, RequestedCoefficientsRetainLocalValues)
{
    basis.nrxx = 1;
    basis.npw = 1;
    basis.nmaxgr = 1;
    charge.nrxx = 1;
    const auto kernel = build("lda", 2, LR::KernelXC::BothSpins);
    ASSERT_EQ(kernel->v3rho3.size(), 4);
    const double expected_singlet = kernel->v3rho3[0] + 2.0 * kernel->v3rho3[1] + kernel->v3rho3[2];
    const double expected_triplet = kernel->v3rho3[0] - 2.0 * kernel->v3rho3[1] + kernel->v3rho3[2];
    ASSERT_EQ(kernel->gxc(false).a_s2.size(), 1);
    ASSERT_EQ(kernel->gxc(true).a_s2.size(), 1);
    EXPECT_DOUBLE_EQ(kernel->gxc(false).a_s2[0], expected_singlet);
    EXPECT_DOUBLE_EQ(kernel->gxc(true).a_s2[0], expected_triplet);
}
