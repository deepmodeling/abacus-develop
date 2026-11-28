#include "../potentials/pot_grad_xc.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_output/cube_io.h"
#include <cmath>
#include <gtest/gtest.h>

// Isolate the real LibXC contraction. The three-point periodic difference and
// divergence below are negative adjoints, so an independently differenced LR
// quadratic form tests the complete local-plus-divergence potential.
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
    const std::complex<double> *in,
    ModuleBase::Vector3<double> *out,
    const ModulePW::PW_Basis *basis,
    const double tpiba)
{
    const ModuleBase::Vector3<double> direction(0.1, 0.2, 0.3);
    for (int point = 0; point < basis->nrxx; ++point)
    {
        const int next = (point + 1) % basis->nrxx;
        const int previous = (point + basis->nrxx - 1) % basis->nrxx;
        const double difference = (in[next].real() - in[previous].real()) * tpiba;
        out[point] = direction * difference;
    }
}
void XC_Functional::grad_dot(
    const ModuleBase::Vector3<double> *in,
    double *out,
    const ModulePW::PW_Basis *basis,
    const double tpiba)
{
    const ModuleBase::Vector3<double> direction(0.1, 0.2, 0.3);
    for (int point = 0; point < basis->nrxx; ++point)
    {
        const int next = (point + 1) % basis->nrxx;
        const int previous = (point + basis->nrxx - 1) % basis->nrxx;
        out[point] = ((in[next] - in[previous]) * direction) * tpiba;
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
    throw std::logic_error("The isolated gradient test does not read coefficient files");
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

/** @brief Compare unrestricted GXC with a separate derivative of the LR quadratic form. */
class PotGradXCTest : public testing::Test
{
  protected:
    ModulePW::PW_Basis basis;
    UnitCell cell;
    Charge charge;
    Parallel_Grid grid;
    std::vector<double> up;
    std::vector<double> down;
    double *density[2];

    void SetUp() override
    {
        basis.nrxx = 3;
        basis.npw = 3;
        basis.nmaxgr = 3;
        cell.omega = 3.0;
        cell.tpiba = 1.0;
        charge.nrxx = 3;
        up = {0.3, 0.3001, 0.2999};
        down = up;
        density[0] = up.data();
        density[1] = down.data();
        charge.rho = density;
    }

    /** @return Whether this installed HSE implementation advertises KXC. */
    bool has_hse_kxc()
    {
        XC_Functional::set_xc_type("hse");
        const double alpha = XC_Functional::get_hybrid_alpha();
        const double omega = XC_Functional::get_hse_omega();
        auto functions = XC_Functional_Libxc::init_func(
            XC_Functional::get_func_id(),
            XC_POLARIZED,
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

    /** @return Borrowed-grid kernel with explicit shell layout and derivative request. */
    std::shared_ptr<const LR::KernelXC> kernel(const std::string &name, const bool open, const int gxc)
    {
        const std::vector<std::string> initialization = {"default"};
        return std::make_shared<const LR::KernelXC>(
            basis,
            cell,
            charge,
            grid,
            2,
            name,
            initialization,
            open,
            gxc);
    }

    /** @return Discrete gradient with the same linear operator as the isolated FFT mock. */
    std::vector<ModuleBase::Vector3<double>> gradient(const std::vector<double> &values) const
    {
        std::vector<std::complex<double>> reciprocal(values.begin(), values.end());
        std::vector<ModuleBase::Vector3<double>> result(3);
        XC_Functional::grad_rho(
            reciprocal.data(),
            result.data(),
            &basis,
            cell.tpiba);
        return result;
    }

    /**
     * @brief Independently evaluate the second-order LR quadratic form, without KXC.
     * @param[in] name Functional for the evaluated LR kernel.
     * @param[in] transition Two real transition-density channels on the three-point grid.
     * @return Unweighted quadratic form in Hartree conventions.
     */
    double quadratic(const std::string &name, const std::vector<double> transition[2])
    {
        const auto second = kernel(name, true, LR::KernelXC::NoGxc);
        const auto first_gradient = gradient(transition[0]);
        const auto second_gradient = gradient(transition[1]);
        double result = 0.0;
        for (int point = 0; point < 3; ++point)
        {
            const double u[2] = {transition[0][point], transition[1][point]};
            const auto &gu = second->drho_gs[0][point];
            const auto &gd = second->drho_gs[1][point];
            const auto &vu = first_gradient[point];
            const auto &vd = second_gradient[point];
            const double linear[3] = {2.0 * (gu * vu), gu * vd + gd * vu, 2.0 * (gd * vd)};
            const double quadratic_gradient[3] = {2.0 * (vu * vu), 2.0 * (vu * vd), 2.0 * (vd * vd)};
            for (int s = 0; s < 2; ++s)
            {
                for (int t = 0; t < 2; ++t)
                {
                    result += second->v2rho2[3 * point + LR::libxc_idx::r2(s, t)] * u[s] * u[t];
                }
                for (int a = 0; a < 3; ++a)
                {
                    result += 2.0 * second->v2rhosigma[6 * point + 3 * s + a] * u[s] * linear[a];
                }
            }
            for (int a = 0; a < 3; ++a)
            {
                result += second->vsigma[3 * point + a] * quadratic_gradient[a];
                for (int b = 0; b < 3; ++b)
                {
                    result += second->v2sigma2[6 * point + LR::libxc_idx::p2[a][b]] * linear[a] * linear[b];
                }
            }
        }
        return result;
    }
};

TEST_F(PotGradXCTest, UnrestrictedReducesToClosedShellWithClippedHse)
{
    if (!has_hse_kxc())
    {
        GTEST_SKIP() << "Installed LibXC does not advertise HSE KXC";
    }
    std::vector<double> transition = {0.8, -0.4, 0.7};
    for (const std::string &name : {"pbe", "hse"})
    {
        const auto closed = kernel(name, false, LR::KernelXC::BothSpins);
        const auto open = kernel(name, true, LR::KernelXC::BothSpins);
        for (const bool triplet : {false, true})
        {
            std::vector<double> second = transition;
            if (triplet)
            {
                for (double &value : second)
                {
                    value = -value;
                }
            }
            double *channels[2] = {transition.data(), second.data()};
            LR::PotGradXCLR reference(
                *closed,
                basis,
                cell,
                3,
                2,
                triplet);
            LR::PotGradXCLR candidate(
                *open,
                basis,
                cell,
                3,
                2,
                false);
            ModuleBase::matrix expected(1, 3);
            reference.cal_v_eff(channels, cell, expected);
            // Closed-shell coefficients pin the free derivative index to up.
            // Compare that same contraction strictly. Both free-spin derivatives
            // are independently checked by the quadratic-form FD test below.
            const int tau = 0;
            ModuleBase::matrix actual(1, 3);
            candidate.cal_v_eff_openshell(
                channels,
                cell,
                actual,
                tau);
            for (int point = 0; point < 3; ++point)
            {
                EXPECT_NEAR(actual(0, point), expected(0, point), 1e-11)
                    << name << " tau=" << tau << " triplet=" << triplet;
            }
        }
    }
}

TEST_F(PotGradXCTest, UnrestrictedHseMatchesIndependentQuadraticFiniteDifference)
{
    if (!has_hse_kxc())
    {
        GTEST_SKIP() << "Installed LibXC does not advertise HSE KXC";
    }
    const std::vector<double> transition[2] = {{0.8, -0.4, 0.7}, {-0.5, 0.6, 0.3}};
    const double perturbation[3] = {0.4, -0.3, 0.2};
    for (const double orientation : {1.0, -1.0})
    {
        up = {0.3, 0.4, 0.2};
        down = {0.18, 0.18 + orientation * 0.005, 0.18 - orientation * 0.005};
        density[0] = up.data();
        density[1] = down.data();
        cell.tpiba = 0.1;
        const auto coefficients = kernel("hse", true, LR::KernelXC::BothSpins);
        const double *channels[2] = {transition[0].data(), transition[1].data()};
        LR::PotGradXCLR potential(
            *coefficients,
            basis,
            cell,
            3,
            2,
            false);
        for (int tau = 0; tau < 2; ++tau)
        {
            ModuleBase::matrix derivative(1, 3);
            potential.cal_v_eff_openshell(
                channels,
                cell,
                derivative,
                tau);
            double analytic = 0.0;
            for (int point = 0; point < 3; ++point)
            {
                analytic += derivative(0, point) * perturbation[point] / ModuleBase::e2;
            }
            std::vector<double> &varied = tau == 0 ? up : down;
            const std::vector<double> original = varied;
            for (const double h : {1e-4, 5e-5})
            {
                for (int point = 0; point < 3; ++point)
                {
                    varied[point] = original[point] + h * perturbation[point];
                }
                const double plus = quadratic("hse", transition);
                for (int point = 0; point < 3; ++point)
                {
                    varied[point] = original[point] - h * perturbation[point];
                }
                const double minus = quadratic("hse", transition);
                for (int point = 0; point < 3; ++point)
                {
                    varied[point] = original[point];
                }
                const double finite_difference = (plus - minus) / (2.0 * h);
                const double tolerance = 2e-8 + 2e-8 * std::abs(finite_difference);
                EXPECT_NEAR(analytic, finite_difference, tolerance)
                    << "tau=" << tau << " orientation=" << orientation << " h=" << h;
            }
        }
    }
}

TEST_F(PotGradXCTest, EmptyLocalGridPreservesEmptySpinChannels)
{
    basis.nrxx = 0;
    basis.npw = 0;
    basis.nmaxgr = 0;
    charge.nrxx = 0;
    const auto coefficients = kernel("pbe", true, LR::KernelXC::BothSpins);
    LR::PotGradXCLR potential(
        *coefficients,
        basis,
        cell,
        0,
        2,
        false);
    const double *channels[2] = {nullptr, nullptr};
    for (int tau = 0; tau < 2; ++tau)
    {
        ModuleBase::matrix result(1, 0);
        EXPECT_NO_THROW(potential.cal_v_eff_openshell(
            channels,
            cell,
            result,
            tau));
        EXPECT_EQ(result.nc, 0);
    }
}


TEST_F(PotGradXCTest, EmptyClosedShellGridRetainsRequestedCoefficients)
{
    basis.nrxx = 0;
    basis.npw = 0;
    basis.nmaxgr = 0;
    charge.nrxx = 0;
    const std::vector<std::string> initialization = {"default"};
    for (const std::string name : {"lda", "pbe"})
    {
        for (const int spins : {1, 2})
        {
            const auto coefficients = LR::build_xc_kernel(name, basis, cell, charge,
                grid, spins, false, LR::KernelXC::BothSpins, initialization);
            double* channels[2] = {nullptr, nullptr};
            for (const bool triplet : {false, true})
            {
                LR::PotGradXCLR potential(*coefficients, basis, cell, 0, spins, triplet);
                ModuleBase::matrix result(1, 0);
                EXPECT_NO_THROW(potential.cal_v_eff(channels, cell, result));
                EXPECT_EQ(result.nc, 0);
            }
        }
    }
}
