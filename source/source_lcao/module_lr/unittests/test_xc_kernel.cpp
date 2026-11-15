#include "../potentials/xc_kernel.h"
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
    const double)
{
    for (int point = 0; point < basis->nrxx; ++point)
    {
        out[point] = ModuleBase::Vector3<double>(0.1, 0.2, 0.3);
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
};

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
