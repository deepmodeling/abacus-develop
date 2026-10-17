#include <gtest/gtest.h>
#include "../kernel_pots.h"
#include "../potentials/pot_grad_xc.h"
#include "source_io/module_parameter/input_parameter.h"
#include "../potentials/pot_hxc_lrtd.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
#include "source_io/module_output/cube_io.h"

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

void XC_Functional::grad_rho(const std::complex<double>*, ModuleBase::Vector3<double>* out,
    const ModulePW::PW_Basis*, const double)
{
    out[0] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
}
void XC_Functional::grad_dot(const ModuleBase::Vector3<double>*, double* out,
    const ModulePW::PW_Basis*, const double)
{
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
    const auto lda = build("lda", 2);
    LR::PotHxcLR original(lda, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    const auto pbe = build("pbe", 2);
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
    const auto pbe = build("pbe", 2);
    LR::PotHxcLR original(pbe, basis, cell, 1, 2, LR::PotHxcLR::S2_triplet);
    const auto lda = build("lda", 2);
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
    for (const auto& name : {"rpa", "hf"})
    {
        const auto kernel = build(name, 2);
        EXPECT_EQ(kernel->nspin(), 2);
        EXPECT_EQ(kernel->functional_type(), LR::PotHxcLR::None);
        EXPECT_EQ(kernel->name(), name);
    }
    init = {"file", "uu", "ud", "dd"};
    const auto file = build("lda", 2);
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
    const auto kernel = build("lda", 1);
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
        const auto prepared = LR::prepare_lr_potentials(input, input.dft_functional,
            basis, cell, charge, grid, nspin, false);
        ASSERT_EQ(prepared.pot.size(), nspin);
        for (const auto& potential : prepared.pot)
        {
            EXPECT_EQ(&potential->xc_kernel_components(), prepared.kernel.get());
        }
        const auto gs = LR::prepare_gs_potential(input, input.dft_functional,
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
    const auto prepared = LR::prepare_lr_potentials(input, lr_functional,
        basis, cell, charge, grid, 2, false);
    const auto gs = LR::prepare_gs_potential(input, lr_functional,
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
    const auto prepared = LR::prepare_lr_potentials(input, input.dft_functional,
        basis, cell, charge, grid, 2, true);
    ASSERT_EQ(prepared.pot.size(), 2);
    EXPECT_EQ(&prepared.pot[0]->xc_kernel_components(), &prepared.pot[1]->xc_kernel_components());
    const auto gs = LR::prepare_gs_potential(input, input.dft_functional,
        basis, cell, charge, grid, 2, true, prepared);
    EXPECT_EQ(&gs->xc_kernel_components(), prepared.kernel.get());
}

TEST_F(PotHxcLRTest, GradientPotentialUsesExplicitKernelSpin)
{
    const auto polarized = build("lda", 2);
    const auto unpolarized = build("lda", 1);
    LR::PotGradXCLR gradient(*polarized, basis, cell, 1, 2, false);
    EXPECT_EQ(gradient.nspin, 2);
    LR::PotGradXCLR single(*unpolarized, basis, cell, 1, 1, false);
    EXPECT_EQ(single.nspin, 1);
    EXPECT_THROW(LR::PotGradXCLR(*polarized, basis, cell, 1, 1, false), std::invalid_argument);
}
