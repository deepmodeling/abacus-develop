#include <gtest/gtest.h>
#include <limits>
#include "source_base/constants.h"
#include "../kernel_pots.h"
#include "../potentials/pot_grad_xc.h"
#include "source_io/module_parameter/input_parameter.h"
#include "../potentials/pot_hxc_lrtd.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
#include "source_io/module_output/cube_io.h"

namespace
{
ModuleBase::Vector3<double> last_gga_divergence;
bool asymmetric_test_gradients = false;
}

// Keep the real kernel/LibXC and potential code. Mock geometry ownership, FFTs,
// external coefficient-file input, and grid masks: exercise binding and metadata.
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
}

namespace ModulePW
{
PW_Basis::PW_Basis() {}
PW_Basis::~PW_Basis() {}
void PW_Basis::initgrids(double, ModuleBase::Matrix3, double) {}
void PW_Basis::initgrids(double, ModuleBase::Matrix3, int, int, int) {}
void PW_Basis::distribute_r() {}
template <typename T>
void PW_Basis::real2recip(const T* in, std::complex<T>* out, const bool, const T) const
{
    out[0] = in[0];
}
template void PW_Basis::real2recip<double>(const double*, std::complex<double>*, bool, double) const;
template <typename T>
void PW_Basis::recip2real(const std::complex<T>* in, std::complex<T>* out, const bool, const T) const
{
    out[0] = in[0];
}
template void PW_Basis::recip2real<double>(const std::complex<double>*, std::complex<double>*, bool, double) const;
}

void XC_Functional::grad_rho(const std::complex<double>* in, ModuleBase::Vector3<double>* out,
    const ModulePW::PW_Basis*, const double)
{
    out[0] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
    if (asymmetric_test_gradients && in[0].real() < 0.25)
    {
        out[0] = (in[0].real() > 0.15)
            ? ModuleBase::Vector3<double>(-0.3, 0.1, 0.05)
            : ModuleBase::Vector3<double>(0.2, -0.15, 0.4);
    }
}
void XC_Functional::grad_dot(const ModuleBase::Vector3<double>* in, double* out,
    const ModulePW::PW_Basis*, const double)
{
    last_gga_divergence = in[0];
    out[0] = 0.;
}
namespace ModuleIO
{
bool read_vdata_palgrid(const Parallel_Grid&, const int, std::ostream&,
    const std::string&, double* const data, const int)
{
    data[0] = 0.25;
    return true;
}
}

namespace XC_Functional_Libxc
{
std::vector<double> cal_sgn(const double, const double, const xc_func_type&,
    const int nspin, const std::size_t nrxx, const std::vector<double>&, const std::vector<double>&)
{
    return std::vector<double>(nrxx * nspin, 1.0);
}
}

class PotHxcLRTest : public testing::Test
{
protected:
    ModulePW::PW_Basis basis;
    UnitCell cell;
    Charge charge;
    Parallel_Grid grid;
    double up = 0.3;
    double down = 0.2;
    double* density[2] = {&up, &down};
    std::vector<std::string> init = {"default"};
    void SetUp() override
    {
        asymmetric_test_gradients = false;
        basis.nrxx = 1;
        basis.npw = 1;
        basis.nmaxgr = 1;
        cell.omega = 1.;
        cell.tpiba = 1.;
        charge.nrxx = 1;
        charge.rho = density;
    }
    std::shared_ptr<const LR::KernelXC> build(const std::string& name, const int nspin)
    {
        return LR::build_xc_kernel(name, basis, cell, charge, grid, nspin, false,
            LR::KernelXC::NoGxc, init);
    }
};

TEST_F(PotHxcLRTest, RebindLdaAfterPbeKeepsMetadataAndSharesCoefficients)
{
    const std::shared_ptr<const LR::KernelXC> lda = build("lda", 2);
    LR::PotHxcLR original(lda, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    const std::shared_ptr<const LR::KernelXC> pbe = build("pbe", 2);
    ASSERT_EQ(XC_Functional::get_func_type(), LR::PotHxcLR::GGA);
    LR::PotHxcLR rebound(lda, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    EXPECT_EQ(lda->functional_type(), LR::PotHxcLR::LDA);
    EXPECT_EQ(pbe->functional_type(), LR::PotHxcLR::GGA);
    EXPECT_EQ(&original.xc_kernel_components(), &rebound.xc_kernel_components());
    ModuleBase::matrix first(1, 1);
    ModuleBase::matrix second(1, 1);
    original.cal_v_eff(density, cell, first);
    rebound.cal_v_eff(density, cell, second);
    EXPECT_DOUBLE_EQ(first(0, 0), second(0, 0));
    EXPECT_NE(first(0, 0), 0.);
    EXPECT_EQ(XC_Functional::get_func_type(), LR::PotHxcLR::GGA);
}

TEST_F(PotHxcLRTest, RebindPbeAfterLdaKeepsGradientPath)
{
    const std::shared_ptr<const LR::KernelXC> pbe = build("pbe", 2);
    LR::PotHxcLR original(pbe, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    const std::shared_ptr<const LR::KernelXC> lda = build("lda", 2);
    LR::PotHxcLR rebound(pbe, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    ModuleBase::matrix first(1, 1);
    ModuleBase::matrix second(1, 1);
    original.cal_v_eff(density, cell, first);
    rebound.cal_v_eff(density, cell, second);
    EXPECT_DOUBLE_EQ(first(0, 0), second(0, 0));
    EXPECT_NE(first(0, 0), 0.);
    EXPECT_EQ(XC_Functional::get_func_type(), LR::PotHxcLR::LDA);
}

TEST_F(PotHxcLRTest, NonlocalAndFileKernelsKeepTheirOwnSpinMetadata)
{
    XC_Functional::set_xc_type("pbe");
    for (const char* const& name : {"rpa", "hf"})
    {
        const std::shared_ptr<const LR::KernelXC> kernel = build(name, 2);
        EXPECT_EQ(kernel->nspin(), 2);
        EXPECT_EQ(kernel->functional_type(), LR::PotHxcLR::None);
        EXPECT_EQ(kernel->name(), name);
    }
    init = {"file", "uu", "ud", "dd"};
    const std::shared_ptr<const LR::KernelXC> file = build("lda", 2);
    EXPECT_EQ(file->nspin(), 2);
    EXPECT_EQ(file->functional_type(), LR::PotHxcLR::LDA);
    EXPECT_EQ(file->v2rho2.size(), 3);
    LR::PotHxcLR potential(file, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    EXPECT_EQ(&potential.xc_kernel_components(), file.get());
}

TEST_F(PotHxcLRTest, RejectMissingOrMismatchedKernel)
{
    std::shared_ptr<const LR::KernelXC> missing;
    EXPECT_THROW(LR::PotHxcLR(missing, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet), std::invalid_argument);
    const std::shared_ptr<const LR::KernelXC> kernel = build("lda", 1);
    EXPECT_EQ(kernel->nspin(), 1);
    EXPECT_THROW(LR::PotHxcLR(kernel, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet), std::invalid_argument);
}

TEST_F(PotHxcLRTest, PotentialFactorySharesSameFunctionalAcrossSpinViews)
{
    Input_para input;
    input.cal_force = true;
    input.dft_functional = "lda";
    input.lr_init_xc_kernel = init;
    for (const int nspin : {1, 2})
    {
        const LR::KernelPotentials prepared = LR::prepare_lr_potentials(input, input.dft_functional,
            basis, cell, charge, grid, nspin, false);
        ASSERT_EQ(prepared.pot.size(), nspin);
        for (const std::shared_ptr<LR::PotHxcLR>& potential : prepared.pot)
        {
            EXPECT_EQ(&potential->xc_kernel_components(), prepared.kernel.get());
        }
        const std::shared_ptr<LR::PotHxcLR> gs = LR::prepare_gs_potential(input, input.dft_functional,
            basis, cell, charge, grid, nspin, false, prepared);
        EXPECT_EQ(&gs->xc_kernel_components(), prepared.kernel.get());
    }
}

TEST_F(PotHxcLRTest, PotentialFactorySeparatesMixedFunctionals)
{
    Input_para input;
    input.cal_force = true;
    input.dft_functional = "pbe";
    input.lr_init_xc_kernel = init;
    const std::string lr_functional = "lda";
    const LR::KernelPotentials prepared = LR::prepare_lr_potentials(input, lr_functional,
        basis, cell, charge, grid, 2, false);
    const std::shared_ptr<LR::PotHxcLR> gs = LR::prepare_gs_potential(input, lr_functional,
        basis, cell, charge, grid, 2, false, prepared);
    EXPECT_NE(&gs->xc_kernel_components(), prepared.kernel.get());
    EXPECT_EQ(gs->xc_kernel_components().name(), "pbe");
    EXPECT_EQ(prepared.kernel->name(), "lda");
}

TEST_F(PotHxcLRTest, PotentialFactoryRetainsOpenShellKernelSharing)
{
    Input_para input;
    input.cal_force = true;
    input.dft_functional = "lda";
    input.lr_init_xc_kernel = init;
    const LR::KernelPotentials prepared = LR::prepare_lr_potentials(input, input.dft_functional,
        basis, cell, charge, grid, 2, true);
    ASSERT_EQ(prepared.pot.size(), 2);
    EXPECT_EQ(&prepared.pot[0]->xc_kernel_components(), &prepared.pot[1]->xc_kernel_components());
    const std::shared_ptr<LR::PotHxcLR> gs = LR::prepare_gs_potential(input, input.dft_functional,
        basis, cell, charge, grid, 2, true, prepared);
    EXPECT_EQ(&gs->xc_kernel_components(), prepared.kernel.get());
}

TEST_F(PotHxcLRTest, GradientPotentialUsesExplicitKernelSpin)
{
    const std::shared_ptr<const LR::KernelXC> polarized = build("lda", 2);
    const std::shared_ptr<const LR::KernelXC> unpolarized = build("lda", 1);
    LR::PotGradXCLR gradient(*polarized, basis, cell, 1, 2, false);
    EXPECT_EQ(gradient.nspin, 2);
    LR::PotGradXCLR single(*unpolarized, basis, cell, 1, 1, false);
    EXPECT_EQ(single.nspin, 1);
    EXPECT_THROW(LR::PotGradXCLR(*polarized, basis, cell, 1, 1, false), std::invalid_argument);
}

TEST_F(PotHxcLRTest, ScopedGsTypeRestoresLrSelection)
{
    const XC_Functional::ScopedType baseline("pbe");
    const std::vector<int> ids = XC_Functional::get_func_id();
    const int type = XC_Functional::get_func_type();
    const std::string description = XC_Functional::output_info();
    {
        const XC_Functional::ScopedType gs("hf");
        EXPECT_TRUE(XC_Functional::get_func_id().empty());
        EXPECT_EQ(XC_Functional::get_func_type(), 4);
        EXPECT_FALSE(XC_Functional::get_ked_flag());
    }
    EXPECT_EQ(XC_Functional::get_func_id(), ids);
    EXPECT_EQ(XC_Functional::get_func_type(), type);
    EXPECT_EQ(XC_Functional::output_info(), description);
}

TEST_F(PotHxcLRTest, ScopedGsTypeRestoresOnException)
{
    const XC_Functional::ScopedType baseline("pbe");
    const std::vector<int> ids = XC_Functional::get_func_id();
    XC_Functional::set_hybrid_alpha(0.37);
    XC_Functional::set_hse_omega(0.19);
    XCFunctionalParameters parameters;
    parameters.xc_temperature = 0.012;
    parameters.exx_fock_alpha = {"0.37"};
    parameters.exx_erfc_alpha = {"0.2"};
    parameters.xc_exch_ext = {1.0, 2.0};
    parameters.xc_corr_ext = {3.0, 4.0};
    XC_Functional::set_runtime_parameters(parameters);
    try
    {
        const XC_Functional::ScopedType gs("hf");
        XC_Functional::set_hybrid_alpha(0.8);
        XC_Functional::set_hse_omega(0.4);
        const XCFunctionalParameters changed;
        XC_Functional::set_runtime_parameters(changed);
        throw std::runtime_error("potential evaluation failed");
    }
    catch (const std::runtime_error&)
    {
    }
    EXPECT_EQ(XC_Functional::get_func_id(), ids);
    EXPECT_EQ(XC_Functional::get_func_type(), 2);
    EXPECT_FALSE(XC_Functional::get_ked_flag());
    EXPECT_FALSE(XC_Functional::get_need_laplacian());
    EXPECT_DOUBLE_EQ(XC_Functional::get_hybrid_alpha(), 0.37);
    EXPECT_DOUBLE_EQ(XC_Functional::get_hse_omega(), 0.19);
    const XCFunctionalParameters& restored = XC_Functional::get_runtime_parameters();
    EXPECT_DOUBLE_EQ(restored.xc_temperature, parameters.xc_temperature);
    EXPECT_EQ(restored.exx_fock_alpha, parameters.exx_fock_alpha);
    EXPECT_EQ(restored.exx_erfc_alpha, parameters.exx_erfc_alpha);
    EXPECT_EQ(restored.xc_exch_ext, parameters.xc_exch_ext);
    EXPECT_EQ(restored.xc_corr_ext, parameters.xc_corr_ext);
}

TEST_F(PotHxcLRTest, ScopedGsTypesNestWithoutLosingEmptyHfSelection)
{
    const XC_Functional::ScopedType baseline("pbe");
    const std::vector<int> ids = XC_Functional::get_func_id();
    {
        const XC_Functional::ScopedType gs("hf");
        {
            const XC_Functional::ScopedType local("lda");
            EXPECT_EQ(XC_Functional::get_func_type(), 1);
            EXPECT_EQ(XC_Functional::get_func_id().size(), 2);
        }
        EXPECT_TRUE(XC_Functional::get_func_id().empty());
        EXPECT_EQ(XC_Functional::get_func_type(), 4);
    }
    EXPECT_EQ(XC_Functional::get_func_id(), ids);
    EXPECT_EQ(XC_Functional::get_func_type(), 2);
}

TEST_F(PotHxcLRTest, ScopedGsTypeRestoresLibxcBackend)
{
    const XC_Functional::ScopedType baseline("XC_LDA_X+XC_LDA_C_PZ");
    const std::vector<int> ids = XC_Functional::get_func_id();
    const std::string description = XC_Functional::output_info();
    ASSERT_NE(description.find("Libxc"), std::string::npos);
    {
        const XC_Functional::ScopedType builtin("pbe");
        EXPECT_EQ(XC_Functional::output_info().find("Libxc"), std::string::npos);
    }
    EXPECT_EQ(XC_Functional::get_func_id(), ids);
    EXPECT_EQ(XC_Functional::output_info(), description);
}

TEST_F(PotHxcLRTest, ScalarTripletMatchesPolarizedClosedShellResponse)
{
    Input_para input;
    input.cal_force = true;
    input.lr_target_spin = "triplet";
    input.lr_init_xc_kernel = init;
    for (const std::string functional : {"lda", "pbe", "hf"})
    {
        input.dft_functional = functional;
        charge.nspin = 1;
        up = 0.6;
        const auto scalar = LR::prepare_lr_potentials(input, functional,
            basis, cell, charge, grid, 1, false);
        ASSERT_EQ(scalar.pot.size(), 1);
        EXPECT_EQ(scalar.kernel->nspin(), 2);
        const auto gs = LR::prepare_gs_potential(input, functional,
            basis, cell, charge, grid, 1, false, scalar);
        EXPECT_EQ(gs->xc_kernel_components().nspin(), 1);
        EXPECT_NE(&gs->xc_kernel_components(), scalar.kernel.get());
        charge.nspin = 2;
        up = 0.3;
        down = 0.3;
        const auto polarized = LR::prepare_lr_potentials(input, functional,
            basis, cell, charge, grid, 2, false);
        EXPECT_EQ(scalar.kernel->v2rho2, polarized.kernel->v2rho2);
        ModuleBase::matrix actual(1, 1);
        ModuleBase::matrix expected(1, 1);
        scalar.pot[0]->cal_v_eff(density, cell, actual);
        polarized.pot[1]->cal_v_eff(density, cell, expected);
        EXPECT_DOUBLE_EQ(actual(0, 0), expected(0, 0));
        if (functional == "hf")
        {
            EXPECT_DOUBLE_EQ(actual(0, 0), 0.0);
            continue;
        }
        // Equal spin densities do not imply f_uu = f_ud: the triplet XC response survives.
        EXPECT_NE(actual(0, 0), 0.0);
        LR::PotGradXCLR scalar_grad(*scalar.kernel, basis, cell, 1, 2, true);
        LR::PotGradXCLR polarized_grad(*polarized.kernel, basis, cell, 1, 2, true);
        actual.zero_out();
        expected.zero_out();
        scalar_grad.cal_v_eff(density, cell, actual);
        polarized_grad.cal_v_eff(density, cell, expected);
        EXPECT_DOUBLE_EQ(actual(0, 0), expected(0, 0));
    }
}

TEST_F(PotHxcLRTest, ScalarTripletRequiresSpinResolvedKernelFiles)
{
    Input_para input;
    input.lr_target_spin = "triplet";
    input.lr_init_xc_kernel = {"file", "scalar"};
    charge.nspin = 1;
    EXPECT_THROW(LR::prepare_lr_potentials(input, "lda", basis, cell, charge, grid, 1, false),
        std::invalid_argument);
    input.lr_init_xc_kernel = {"file", "uu", "ud", "dd"};
    const auto prepared = LR::prepare_lr_potentials(input, "lda", basis, cell, charge, grid, 1, false);
    EXPECT_EQ(prepared.kernel->v2rho2.size(), 3);
}

TEST_F(PotHxcLRTest, ClosedShellScalarReuseMatchesIndependentGxcContractions)
{
    down = up;
    for (const int nspin : {1, 2})
    {
        const auto kernel = LR::build_xc_kernel("pbe", basis, cell, charge, grid,
            nspin, false, LR::KernelXC::BothSpins, init);
        const int spin_cases = (nspin == 1) ? 1 : 2;
        for (int spin = 0; spin < spin_cases; ++spin)
        {
            const bool triplet = (spin == 1);
            const auto& scalar = kernel->closed_shell_gga(triplet);
            ASSERT_EQ(scalar.rho_sigma.size(), 1u);
            ASSERT_EQ(scalar.sigma_sigma.size(), 1u);
            const auto& g = kernel->gxc(triplet);
            const double eta[2] = {1., triplet ? -1. : 1.};
            const double mu[3] = {1., triplet ? 0. : 1., triplet ? -1. : 1.};
            const double weight[3] = {2., triplet ? -1. : 1., 0.};
            double cs = 0.;
            double ct = 0.;
            if (nspin == 1)
            {
                EXPECT_DOUBLE_EQ(scalar.rho_sigma[0], 2. * kernel->v2rhosigma[0]);
                EXPECT_DOUBLE_EQ(scalar.sigma_sigma[0], 4. * kernel->v2sigma2[0]);
                cs = 16. * kernel->v2rhosigma[0];
                ct = 32. * kernel->v2sigma2[0];
            }
            else
            {
                for (int a = 0; a < 3; ++a)
                {
                    for (int b = 0; b < 2; ++b)
                    {
                        cs += 2. * weight[a] * eta[b] * kernel->v2rhosigma[3 * b + a];
                    }
                    for (int b = 0; b < 3; ++b)
                    {
                        ct += 4. * weight[a] * mu[b] * kernel->v2sigma2[LR::libxc_idx::p2[a][b]];
                    }
                }
                EXPECT_NEAR(2. * scalar.rho_sigma[0], cs, 1e-11);
                EXPECT_NEAR(2. * scalar.sigma_sigma[0], ct, 1e-11);
            }
            LR::PotGradXCLR potential(*kernel, basis, cell, 1, nspin, triplet);
            ModuleBase::matrix result(1, 1);
            potential.cal_v_eff(density, cell, result);
            const auto& dg = kernel->drho_gs[0][0];
            const ModuleBase::Vector3<double> dx(0.1, 0.2, 0.3);
            const double t = dg * dx;
            const double q = dx * dx;
            const double e = g.e_s2[0] * up * up + g.e_st[0] * up * t
                + g.e_t2[0] * t * t + g.e_q[0] * q;
            const auto expected = -(dg * e + dx * (cs * up + ct * t));
            EXPECT_NEAR(last_gga_divergence.x, expected.x, 1e-11);
            EXPECT_NEAR(last_gga_divergence.y, expected.y, 1e-11);
            EXPECT_NEAR(last_gga_divergence.z, expected.z, 1e-11);
        }
    }
}

TEST_F(PotHxcLRTest, HseGroundResponseDoesNotReuseLrSigmaFloor)
{
    Input_para input;
    input.dft_functional = "hse";
    input.lr_init_xc_kernel = init;
    for (const int spins : {1, 2})
    {
        const auto prepared = LR::prepare_lr_potentials(input, input.dft_functional,
            basis, cell, charge, grid, spins, false);
        const auto gs = LR::prepare_gs_potential(input, input.dft_functional,
            basis, cell, charge, grid, spins, false, prepared);
        EXPECT_NE(&gs->xc_kernel_components(), prepared.kernel.get());
        EXPECT_EQ(gs->xc_kernel_components().nspin(), spins);
        EXPECT_TRUE(gs->xc_kernel_components().v3rho3.empty());
    }
}

TEST_F(PotHxcLRTest, UnrestrictedScalarContractionPreservesAllSpinBlocks)
{
    asymmetric_test_gradients = true;
    basis.ig_gge0 = 0;
    const auto kernel = LR::build_xc_kernel("pbe", basis, cell, charge, grid,
        2, true, LR::KernelXC::BothSpins, init);
    LR::PotHxcLR potential(kernel, basis, cell, 1, 2, LR::PotHxcLR::S2_updown);
    const auto& gu = kernel->drho_gs[0][0];
    const auto& gd = kernel->drho_gs[1][0];
    const auto gradient_cross = ModuleBase::cross(gu, gd);
    EXPECT_GT(gradient_cross.norm(), 0.01);
    const ModuleBase::Vector3<double> zero(0., 0., 0.);
    const ModuleBase::Vector3<double> sigma_gradient[2][3] = {
        {gu * 2., gd, zero}, {zero, gu, gd * 2.}};
    const ModuleBase::Vector3<double> perturbation_gradient(0.2, -0.15, 0.4);
    double perturbation = 0.07;
    double* response[1] = {&perturbation};
    for (int output = 0; output < 2; ++output)
    {
        for (int input = 0; input < 2; ++input)
        {
            ModuleBase::matrix result(1, 1);
            const std::vector<int> block = {output, input};
            potential.cal_v_eff(response, cell, result, block);
            ModuleBase::Vector3<double> mixed_divergence(0., 0., 0.);
            ModuleBase::Vector3<double> mixed_local(0., 0., 0.);
            ModuleBase::Vector3<double> hessian_action(0., 0., 0.);
            for (int a = 0; a < 3; ++a)
            {
                mixed_divergence += sigma_gradient[output][a] * kernel->v2rhosigma[3 * input + a];
                mixed_local += sigma_gradient[input][a] * kernel->v2rhosigma[3 * output + a];
                for (int b = 0; b < 3; ++b)
                {
                    const double coefficient = kernel->v2sigma2[LR::libxc_idx::p2[a][b]];
                    const double contraction = sigma_gradient[input][b] * perturbation_gradient;
                    hessian_action += sigma_gradient[output][a] * (coefficient * contraction);
                }
            }
            const int component = (output == input) ? 2 * output : 1;
            const double multiplicity = (output == input) ? 2. : 1.;
            hessian_action += perturbation_gradient * (multiplicity * kernel->vsigma[component]);
            const auto expected = -(mixed_divergence * perturbation + hessian_action);
            EXPECT_NEAR(last_gga_divergence.x, expected.x, 1e-12);
            EXPECT_NEAR(last_gga_divergence.y, expected.y, 1e-12);
            EXPECT_NEAR(last_gga_divergence.z, expected.z, 1e-12);
            const double local = kernel->v2rho2[output + input] * perturbation
                + mixed_local * perturbation_gradient;
            EXPECT_NEAR(result(0, 0), ModuleBase::e2 * local, 1e-12);
        }
    }
    asymmetric_test_gradients = false;
}

TEST_F(PotHxcLRTest, UnrestrictedBlockDoesNotReadInactiveScalarComponents)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double rho_sigma[6] = {1., 2., 3., nan, nan, nan};
    const double sigma_sigma[6] = {4., 5., nan, 6., nan, nan};
    const double density_response[2] = {0.2, 0.};
    const double sigma_response[3] = {0.1, 0.3, 0.};
    const double first = LR::sigma_potential_derivative(
        rho_sigma, sigma_sigma, 0, density_response, sigma_response, 0, 1);
    const double second = LR::sigma_potential_derivative(
        rho_sigma, sigma_sigma, 1, density_response, sigma_response, 0, 1);
    EXPECT_DOUBLE_EQ(first, 0.2 + 0.4 + 1.5);
    EXPECT_DOUBLE_EQ(second, 0.4 + 0.5 + 1.8);
}
