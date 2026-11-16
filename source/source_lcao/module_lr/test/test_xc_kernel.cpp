#include "../potentials/xc_kernel.h"
#include "../utils/lr_util.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_output/cube_io.h"
#include <gtest/gtest.h>

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
