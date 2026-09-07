#include "../libxc_abacus.h"
#include "../xc_functional.h"
#include "../xc_functional_ncgga_sf.h"
#include "../xc_ncgga_radial.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"

#ifdef __MPI
#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"

#include <mpi.h>
#endif

#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

// This focused target does not link the full elecstate object library. The
// fixture supplies vector-backed charge storage, so only the trivial lifetime
// boundary is needed here.
Charge::Charge()
{
}
Charge::~Charge()
{
}

// This target links the production PW/XC objects but not the full cell object
// library.  The stress entry point only reads UnitCell::tpiba and magnet.lsign_
// in the parent implementation, so provide the same narrow lifetime boundary
// used by the existing XC focused tests.
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}

namespace
{
int test_rank = 0;
int test_size = 1;

double pool_sum(const double local)
{
    double global = local;
    Parallel_Reduce::reduce_pool(global);
    return global;
}

double pool_min(const double local)
{
    double global = local;
    Parallel_Reduce::reduce_min(global);
    return global;
}

double pool_max(const double local)
{
    double global = local;
    Parallel_Reduce::reduce_max(global);
    return global;
}

bool is_pool_root()
{
    return test_rank == 0;
}

std::uint64_t potential_fnv1a64(const ModuleBase::matrix& potential)
{
    static_assert(sizeof(double) == sizeof(std::uint64_t), "the potential hash requires 64-bit doubles");
    std::uint64_t hash = 14695981039346656037ULL;
    for (int channel = 0; channel < potential.nr; ++channel)
    {
        for (int ir = 0; ir < potential.nc; ++ir)
        {
            std::uint64_t bits = 0;
            const double value = potential(channel, ir);
            std::memcpy(&bits, &value, sizeof(bits));
            for (int byte = 0; byte < 8; ++byte)
            {
                hash ^= (bits >> (8 * byte)) & 0xffULL;
                hash *= 1099511628211ULL;
            }
        }
    }
    return hash;
}

class RealPwNcgga : public testing::Test
{
  protected:
    typedef std::tuple<double, double, ModuleBase::matrix> VxcResult;
    typedef std::function<VxcResult()> Evaluator;
    typedef std::function<std::vector<double>()> StressEvaluator;

    struct BranchMargins
    {
        double min_abs_total_density;
        double min_signed_saturation_gap;
        double max_signed_saturation_gap;
        double min_magnitude;
        double max_magnitude;
        double min_eta_distance;
    };

    ModulePW::PW_Basis pw;
    Charge charge;
    std::array<std::vector<double>, 4> density;
    std::array<std::vector<double>, 4> perturbation;
    std::array<double*, 4> density_pointer;
    std::vector<double> core_density;
    std::vector<std::complex<double>> core_density_reciprocal;
    std::array<std::vector<std::complex<double>>, 2> charge_reciprocal;
    std::array<std::complex<double>*, 2> charge_reciprocal_pointer;

    struct ReciprocalMetricState
    {
        std::vector<ModuleBase::Vector3<double>> gcar;
        std::array<std::vector<double>, 4> density;
        std::vector<double> core_density;
        std::vector<std::complex<double>> core_density_reciprocal;
        double omega = 0.0;
    };

    void SetUp() override
    {
#ifdef __MPI
        pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
#endif
        // Keep tpiba away from one so an omitted or duplicated reciprocal-
        // length factor cannot accidentally pass the adjoint test.
        const double lat0 = 7.0;
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        // An odd z dimension gives unequal real-space slabs in the MPI2 test.
        pw.initgrids(lat0, lattice, 24, 10, 9);
        pw.initparameters(false, 80.0, 2, false);
        pw.setuptransform();
        pw.collect_local_pw();

        ASSERT_EQ(pw.nx, 24);
        ASSERT_EQ(pw.ny, 10);
        ASSERT_EQ(pw.nz, 9);
        ASSERT_EQ(pw.nxyz, 2160);
        ASSERT_GT(pw.npwtot, 100);
        ASSERT_FALSE(pw.gamma_only);
        ASSERT_NEAR(pw.tpiba, ModuleBase::TWO_PI / lat0, 1.0e-14);
        ASSERT_GT(std::abs(pw.tpiba - 1.0), 5.0e-2);

        for (int channel = 0; channel < 4; ++channel)
        {
            density[channel].resize(pw.nrxx);
            perturbation[channel].resize(pw.nrxx);
            density_pointer[channel] = density[channel].data();
        }
        core_density.resize(pw.nrxx);
        core_density_reciprocal.resize(pw.npw);
        for (int spin = 0; spin < 2; ++spin)
        {
            charge_reciprocal[spin].resize(pw.npw);
            charge_reciprocal_pointer[spin] = charge_reciprocal[spin].data();
        }

        charge.rhopw = &pw;
        charge.nrxx = pw.nrxx;
        charge.nxyz = pw.nxyz;
        charge.ngmc = pw.npw;
        charge.nspin = 4;
        charge.rho = density_pointer.data();
        charge.rhog = charge_reciprocal_pointer.data();
        charge.rho_core = core_density.data();
        charge.rhog_core = core_density_reciprocal.data();

        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            // PW_Basis stores local real data as
            // ir = iz_local + (iy + ix * ny) * nplane.
            const int ix = ir / (pw.ny * pw.nplane);
            const int iy = (ir / pw.nplane) % pw.ny;
            const int iz = ir % pw.nplane + pw.startz_current;
            const double x = ModuleBase::TWO_PI * static_cast<double>(ix) / pw.nx;
            const double y = ModuleBase::TWO_PI * static_cast<double>(iy) / pw.ny;
            const double z = ModuleBase::TWO_PI * static_cast<double>(iz) / pw.nz;
            const double total_density = 2.2 + 0.18 * std::sin(x + 0.21) + 0.13 * std::cos(4.0 * x - 0.17)
                                         + 0.07 * std::sin(8.0 * x + 0.33) + 0.05 * std::cos(y - 0.26)
                                         + 0.04 * std::sin(z + 0.31);
            const double magnitude = 0.62 + 0.07 * std::cos(2.0 * x + 0.13) + 0.05 * std::sin(5.0 * x - 0.27)
                                     + 0.03 * std::sin(y + 0.19) + 0.02 * std::cos(z - 0.23);
            const double theta = 0.7 + 0.32 * std::sin(3.0 * x + 0.11) + 0.18 * std::cos(7.0 * x - 0.23)
                                 + 0.10 * std::cos(y + 0.17) + 0.07 * std::sin(z - 0.29);
            const double phi = 0.4 + 0.27 * std::cos(4.0 * x + 0.37) - 0.16 * std::sin(6.0 * x + 0.29)
                               + 0.09 * std::sin(y + z + 0.15);
            density[0][ir] = total_density;
            density[1][ir] = magnitude * std::sin(theta) * std::cos(phi);
            density[2][ir] = magnitude * std::sin(theta) * std::sin(phi);
            density[3][ir] = magnitude * std::cos(theta);
            core_density[ir] = 0.15 + 0.03 * std::cos(3.0 * x - 0.14) + 0.02 * std::sin(6.0 * x + 0.25)
                               + 0.01 * std::cos(y - z + 0.18);

            perturbation[0][ir]
                = 0.31 * std::cos(2.0 * x + 0.41) - 0.19 * std::sin(7.0 * x - 0.12) + 0.11 * std::cos(y + z - 0.16);
            perturbation[1][ir]
                = 0.27 * std::sin(x + 0.37) + 0.21 * std::cos(8.0 * x + 0.19) + 0.13 * std::sin(y - 0.24);
            perturbation[2][ir]
                = 0.29 * std::cos(3.0 * x - 0.22) - 0.17 * std::sin(6.0 * x + 0.31) + 0.12 * std::cos(z + 0.28);
            perturbation[3][ir]
                = 0.25 * std::sin(5.0 * x + 0.18) + 0.23 * std::cos(7.0 * x - 0.29) + 0.10 * std::sin(y - z + 0.32);
        }
        pw.real2recip(core_density.data(), core_density_reciprocal.data());
    }

    VxcResult evaluate_builtin(const std::string& functional = "PBE")
    {
        XC_Functional::set_xc_type(functional);
        UnitCell cell;
        cell.omega = pw.omega;
        cell.tpiba = pw.tpiba;
        return XC_Functional::v_xc(pw.nrxx,
                                   &charge,
                                   &cell,
                                   4,
                                   true,
                                   false,
                                   2,
                                   XC_Functional::get_hybrid_alpha(),
                                   XC_Functional::get_hse_omega());
    }

#ifdef __LIBXC
    VxcResult evaluate_libxc(const std::vector<int>& functionals, const std::map<int, double>* scaling_factor)
    {
        return XC_Functional_Libxc::v_xc_libxc(functionals,
                                               pw.nrxx,
                                               pw.omega,
                                               pw.tpiba,
                                               &charge,
                                               4,
                                               true,
                                               false,
                                               2,
                                               scaling_factor,
                                               0.0,
                                               0.0);
    }

    VxcResult evaluate_libxc_gga(const std::map<int, double>* scaling_factor = nullptr)
    {
        const std::vector<int> functionals = {XC_GGA_X_PBE, XC_GGA_C_PBE};
        return evaluate_libxc(functionals, scaling_factor);
    }

    VxcResult evaluate_libxc_lda()
    {
        const std::vector<int> functionals = {XC_LDA_X, XC_LDA_C_PZ};
        return evaluate_libxc(functionals, nullptr);
    }

#endif

    void set_uniform_state(const double total_density, const std::array<double, 3>& magnetization)
    {
        const std::array<double, 4> constant_direction = {{0.17, -0.11, 0.13, 0.09}};
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            density[0][ir] = total_density;
            for (int mu = 0; mu < 3; ++mu)
            {
                density[mu + 1][ir] = magnetization[mu];
            }
            core_density[ir] = 0.0;
            for (int channel = 0; channel < 4; ++channel)
            {
                perturbation[channel][ir] += constant_direction[channel];
            }
        }
        std::fill(core_density_reciprocal.begin(), core_density_reciprocal.end(), std::complex<double>(0.0, 0.0));
    }

    void set_inside_eta_state()
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const int ix = ir / (pw.ny * pw.nplane);
            const int iy = (ir / pw.nplane) % pw.ny;
            const int iz = ir % pw.nplane + pw.startz_current;
            const double x = ModuleBase::TWO_PI * static_cast<double>(ix) / pw.nx;
            const double y = ModuleBase::TWO_PI * static_cast<double>(iy) / pw.ny;
            const double z = ModuleBase::TWO_PI * static_cast<double>(iz) / pw.nz;
            density[0][ir] = 0.030 + 0.002 * std::cos(x - y + 0.2);
            density[1][ir] = 2.8e-4 + 0.6e-4 * std::sin(2.0 * x + 0.1);
            density[2][ir] = -2.1e-4 + 0.5e-4 * std::cos(3.0 * x - z + 0.3);
            density[3][ir] = 1.7e-4 + 0.4e-4 * std::sin(y + z - 0.2);
            core_density[ir] = 0.004 + 0.001 * std::cos(2.0 * x + z - 0.1);
            perturbation[0][ir] += 0.07;
            perturbation[1][ir] += 0.13;
            perturbation[2][ir] -= 0.11;
            perturbation[3][ir] += 0.09;
        }
        pw.real2recip(core_density.data(), core_density_reciprocal.data());
    }

    void set_negative_gga_state()
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const int ix = ir / (pw.ny * pw.nplane);
            const int iy = (ir / pw.nplane) % pw.ny;
            const int iz = ir % pw.nplane + pw.startz_current;
            const double x = ModuleBase::TWO_PI * static_cast<double>(ix) / pw.nx;
            const double y = ModuleBase::TWO_PI * static_cast<double>(iy) / pw.ny;
            const double z = ModuleBase::TWO_PI * static_cast<double>(iz) / pw.nz;
            density[0][ir] = -1.55 - 0.10 * std::cos(x - y + 0.2) - 0.05 * std::sin(3.0 * x + z - 0.1);
            density[1][ir] = 0.30 + 0.05 * std::sin(2.0 * x + 0.1);
            density[2][ir] = -0.24 + 0.04 * std::cos(3.0 * x - z + 0.3);
            density[3][ir] = 0.20 + 0.03 * std::sin(y + z - 0.2);
            core_density[ir] = -0.12 - 0.02 * std::cos(2.0 * x + z - 0.1);
        }
        pw.real2recip(core_density.data(), core_density_reciprocal.data());
    }

    void set_saturated_gga_state()
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const int ix = ir / (pw.ny * pw.nplane);
            const int iy = (ir / pw.nplane) % pw.ny;
            const int iz = ir % pw.nplane + pw.startz_current;
            const double x = ModuleBase::TWO_PI * static_cast<double>(ix) / pw.nx;
            const double y = ModuleBase::TWO_PI * static_cast<double>(iy) / pw.ny;
            const double z = ModuleBase::TWO_PI * static_cast<double>(iz) / pw.nz;
            density[0][ir] = 0.43 + 0.04 * std::cos(x - y + 0.2) + 0.02 * std::sin(3.0 * x + z - 0.1);
            density[1][ir] = 0.64 + 0.06 * std::sin(2.0 * x + 0.1);
            density[2][ir] = 0.34 + 0.05 * std::cos(3.0 * x - z + 0.3);
            density[3][ir] = 0.28 + 0.04 * std::sin(y + z - 0.2);
            core_density[ir] = 0.05 + 0.01 * std::cos(2.0 * x + z - 0.1);
        }
        pw.real2recip(core_density.data(), core_density_reciprocal.data());
    }

    void set_zero_core_density()
    {
        std::fill(core_density.begin(), core_density.end(), 0.0);
        std::fill(core_density_reciprocal.begin(), core_density_reciprocal.end(), std::complex<double>(0.0, 0.0));
    }

    ReciprocalMetricState capture_reciprocal_metric_state() const
    {
        ReciprocalMetricState state;
        state.gcar.assign(pw.gcar, pw.gcar + pw.npw);
        state.density = density;
        state.core_density = core_density;
        state.core_density_reciprocal = core_density_reciprocal;
        state.omega = pw.omega;
        return state;
    }

    void restore_reciprocal_metric_state(const ReciprocalMetricState& state)
    {
        ASSERT_EQ(state.gcar.size(), static_cast<std::size_t>(pw.npw));
        for (int ig = 0; ig < pw.npw; ++ig)
        {
            pw.gcar[ig] = state.gcar[ig];
        }
        for (int channel = 0; channel < 4; ++channel)
        {
            ASSERT_EQ(state.density[channel].size(), density[channel].size());
            std::copy(state.density[channel].begin(), state.density[channel].end(), density[channel].begin());
        }
        ASSERT_EQ(state.core_density.size(), core_density.size());
        std::copy(state.core_density.begin(), state.core_density.end(), core_density.begin());
        ASSERT_EQ(state.core_density_reciprocal.size(), core_density_reciprocal.size());
        std::copy(state.core_density_reciprocal.begin(),
                  state.core_density_reciprocal.end(),
                  core_density_reciprocal.begin());
        pw.omega = state.omega;
    }

    static void add_matrix_element(ModuleBase::Matrix3& matrix, const int row, const int column, const double value)
    {
        ASSERT_GE(row, 0);
        ASSERT_LT(row, 3);
        ASSERT_GE(column, 0);
        ASSERT_LT(column, 3);
        double* element = nullptr;
        if (row == 0 && column == 0)
            element = &matrix.e11;
        if (row == 0 && column == 1)
            element = &matrix.e12;
        if (row == 0 && column == 2)
            element = &matrix.e13;
        if (row == 1 && column == 0)
            element = &matrix.e21;
        if (row == 1 && column == 1)
            element = &matrix.e22;
        if (row == 1 && column == 2)
            element = &matrix.e23;
        if (row == 2 && column == 0)
            element = &matrix.e31;
        if (row == 2 && column == 1)
            element = &matrix.e32;
        if (row == 2 && column == 2)
            element = &matrix.e33;
        ASSERT_NE(element, nullptr);
        *element += value;
    }

    double evaluate_reciprocal_metric_deformation(const ReciprocalMetricState& state,
                                                  const int stress_row,
                                                  const int stress_column,
                                                  const double epsilon,
                                                  const bool homogeneous_density_scaling,
                                                  const Evaluator& evaluator)
    {
        restore_reciprocal_metric_state(state);

        // gcar is stored as a row vector.  Under r' = F r, reciprocal
        // vectors transform as k' = k F^{-1}.  Perturb F_{column,row}; its
        // derivative contracts exactly with the lower-triangle convention
        // sum_r h_row * g_column used by production stress_gga.
        ModuleBase::Matrix3 deformation;
        add_matrix_element(deformation, stress_column, stress_row, epsilon);
        const double determinant = deformation.Det();
        const ModuleBase::Matrix3 inverse = deformation.Inverse();
        for (int ig = 0; ig < pw.npw; ++ig)
        {
            pw.gcar[ig] = state.gcar[ig] * inverse;
        }

        if (homogeneous_density_scaling)
        {
            pw.omega = state.omega * determinant;
            for (int channel = 0; channel < 4; ++channel)
            {
                for (int ir = 0; ir < pw.nrxx; ++ir)
                {
                    density[channel][ir] = state.density[channel][ir] / determinant;
                }
            }
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] = state.core_density[ir] / determinant;
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
        }

        const double energy = std::get<0>(evaluator());
        restore_reciprocal_metric_state(state);
        return energy;
    }

    std::vector<double> evaluate_gradient_stress_dispatch(const std::string& functional)
    {
        XC_Functional::set_xc_type(functional);
        UnitCell cell;
        cell.tpiba = pw.tpiba;
        cell.magnet.lsign_ = false;
        double dummy_energy = 0.0;
        double dummy_vtxc = 0.0;
        ModuleBase::matrix dummy_potential;
        std::vector<double> stress;
        XC_Functional::gradcorr(dummy_energy,
                                dummy_vtxc,
                                dummy_potential,
                                &charge,
                                &pw,
                                &cell,
                                stress,
                                true,
                                4,
                                true,
                                false,
                                2,
                                0.0,
                                0.0);
        EXPECT_EQ(stress.size(), 9U);
        return stress;
    }

    std::vector<double> evaluate_builtin_gradient_stress()
    {
        return evaluate_gradient_stress_dispatch("PBE");
    }

#ifdef __LIBXC
    std::vector<double> evaluate_libxc_pbe_gradient_stress_dispatch()
    {
        return evaluate_gradient_stress_dispatch("GGA_X_PBE+GGA_C_PBE");
    }
#endif

    void expect_gradient_stress_metric_derivative(const std::string& mode,
                                                  const Evaluator& energy_evaluator,
                                                  const StressEvaluator& stress_evaluator,
                                                  const double relative_tolerance,
                                                  const double absolute_floor)
    {
        const ReciprocalMetricState state = capture_reciprocal_metric_state();
        const std::vector<double> local_stress = stress_evaluator();
        const std::array<double, 4> eps_values = {{2.0e-3, 1.0e-3, 5.0e-4, 2.5e-4}};

        for (int row = 0; row < 3; ++row)
        {
            for (int column = 0; column <= row; ++column)
            {
                SCOPED_TRACE(mode + " metric component " + std::to_string(row) + std::to_string(column));
                const int index = row * 3 + column;
                const double analytic = pool_sum(local_stress[index]) / pw.nxyz;
                std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};
                for (int ieps = 0; ieps < 4; ++ieps)
                {
                    const double epsilon = eps_values[ieps];
                    const double energy_plus
                        = evaluate_reciprocal_metric_deformation(state, row, column, epsilon, false, energy_evaluator);
                    const double energy_minus
                        = evaluate_reciprocal_metric_deformation(state, row, column, -epsilon, false, energy_evaluator);
                    const double finite_difference = -(energy_plus - energy_minus) / (2.0 * epsilon * state.omega);
                    errors[ieps] = std::abs(finite_difference - analytic);
                    if (is_pool_root())
                    {
                        const double relative_scale
                            = std::max(1.0e-30, std::max(std::abs(analytic), std::abs(finite_difference)));
                        const double order = ieps == 0 || errors[ieps] == 0.0
                                                 ? 0.0
                                                 : std::log(errors[ieps - 1] / errors[ieps]) / std::log(2.0);
                        std::cout << std::setprecision(17) << "NCGGA_STRESS_METRIC_FD mode=" << mode << " row=" << row
                                  << " column=" << column << " eps=" << epsilon << " analytic=" << analytic
                                  << " finite_difference=" << finite_difference << " absolute_error=" << errors[ieps]
                                  << " relative_error=" << errors[ieps] / relative_scale
                                  << " convergence_order=" << order << '\n';
                    }
                }
                // Stress components can be much smaller than one.  Scaling
                // the tolerance by max(1, |sigma|) would hide a persistent
                // O(10%) slope offset on the regularized radial branch.
                // The absolute floor only covers the observed FFT/energy
                // roundoff plateau; the relative term still resolves every
                // nonzero component in this fixture.
                const double scale = std::max(1.0e-10, std::abs(analytic));
                const double tolerance = relative_tolerance * scale + absolute_floor;
                EXPECT_LE(*std::min_element(errors.begin(), errors.end()), tolerance);
                // Once the absolute error is at the MPI energy-roundoff
                // plateau, individual refinements need not remain monotone.
                // Require the central-difference O(eps^2) contraction only
                // while both adjacent errors are still resolved above that
                // plateau.  The closure check above remains unconditional.
                const double roundoff_plateau = 10.0 * absolute_floor;
                if (errors[0] > roundoff_plateau && errors[1] > roundoff_plateau)
                {
                    EXPECT_LE(errors[1], 0.4 * errors[0] + tolerance);
                }
                if (errors[1] > roundoff_plateau && errors[2] > roundoff_plateau)
                {
                    EXPECT_LE(errors[2], 0.4 * errors[1] + tolerance);
                }
            }
        }
    }

    void expect_builtin_gradient_stress_metric_derivative(const std::string& mode)
    {
        expect_gradient_stress_metric_derivative(
            mode,
            [this]() { return evaluate_builtin(); },
            [this]() { return evaluate_builtin_gradient_stress(); },
            2.0e-5,
            2.0e-13);
    }

    void expect_full_xc_diagonal_stress(const std::string& mode,
                                        const Evaluator& energy_evaluator,
                                        const StressEvaluator& stress_evaluator)
    {
        // The returned vtxc is a valence-only four-channel inner product.
        // Core deformation is accounted separately by stress_cc in the full
        // PW stress path, so this isolated XC diagonal identity requires a
        // zero core density.
        set_zero_core_density();
        ASSERT_EQ(pool_max(*std::max_element(core_density.begin(), core_density.end())), 0.0);
        const ReciprocalMetricState state = capture_reciprocal_metric_state();
        const VxcResult base = energy_evaluator();
        const std::vector<double> local_stress = stress_evaluator();
        const double diagonal_local_term = -(std::get<0>(base) - std::get<1>(base)) / state.omega;
        const std::array<double, 4> eps_values = {{2.0e-3, 1.0e-3, 5.0e-4, 2.5e-4}};

        for (int diagonal = 0; diagonal < 3; ++diagonal)
        {
            SCOPED_TRACE(mode + " full diagonal " + std::to_string(diagonal));
            const double gradient_correction = pool_sum(local_stress[diagonal * 3 + diagonal]) / pw.nxyz;
            const double analytic = diagonal_local_term + gradient_correction;
            std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};
            for (int ieps = 0; ieps < 4; ++ieps)
            {
                const double epsilon = eps_values[ieps];
                const double energy_plus = evaluate_reciprocal_metric_deformation(state,
                                                                                  diagonal,
                                                                                  diagonal,
                                                                                  epsilon,
                                                                                  true,
                                                                                  energy_evaluator);
                const double energy_minus = evaluate_reciprocal_metric_deformation(state,
                                                                                   diagonal,
                                                                                   diagonal,
                                                                                   -epsilon,
                                                                                   true,
                                                                                   energy_evaluator);
                const double finite_difference = -(energy_plus - energy_minus) / (2.0 * epsilon * state.omega);
                errors[ieps] = std::abs(finite_difference - analytic);
                if (is_pool_root())
                {
                    const double relative_scale
                        = std::max(1.0e-30, std::max(std::abs(analytic), std::abs(finite_difference)));
                    const double order = ieps == 0 || errors[ieps] == 0.0
                                             ? 0.0
                                             : std::log(errors[ieps - 1] / errors[ieps]) / std::log(2.0);
                    std::cout << std::setprecision(17) << "NCGGA_STRESS_FULL_FD mode=" << mode
                              << " diagonal=" << diagonal << " eps=" << epsilon
                              << " gradient_correction=" << gradient_correction
                              << " local_diagonal=" << diagonal_local_term << " analytic=" << analytic
                              << " finite_difference=" << finite_difference << " absolute_error=" << errors[ieps]
                              << " relative_error=" << errors[ieps] / relative_scale << " convergence_order=" << order
                              << '\n';
                }
            }
            const double scale = std::max(1.0, std::abs(analytic));
            EXPECT_LE(*std::min_element(errors.begin(), errors.end()), 5.0e-8 * scale);
            EXPECT_LE(errors[1], 0.4 * errors[0] + 1.0e-8 * scale);
            EXPECT_LE(errors[2], 0.4 * errors[1] + 1.0e-8 * scale);
        }
    }

    void expect_builtin_full_xc_diagonal_stress(const std::string& mode)
    {
        expect_full_xc_diagonal_stress(
            mode,
            [this]() { return evaluate_builtin(); },
            [this]() { return evaluate_builtin_gradient_stress(); });
    }

    BranchMargins report_branch_margins(const std::string& mode)
    {
        constexpr double lca_eta = 1.0e-3;
        double local_min_abs_density = std::numeric_limits<double>::max();
        double local_min_gap = std::numeric_limits<double>::max();
        double local_max_gap = -std::numeric_limits<double>::max();
        double local_min_magnitude = std::numeric_limits<double>::max();
        double local_max_magnitude = 0.0;
        double local_min_eta_distance = std::numeric_limits<double>::max();
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const double total = density[0][ir] + core_density[ir];
            const double magnitude = std::sqrt(density[1][ir] * density[1][ir] + density[2][ir] * density[2][ir]
                                               + density[3][ir] * density[3][ir]);
            const ModuleXC::NcggaRadialPoint radial
                = ModuleXC::make_ncgga_radial_point({{density[1][ir], density[2][ir], density[3][ir]}}, lca_eta);
            const double gap = std::abs(total) - radial.value;
            local_min_abs_density = std::min(local_min_abs_density, std::abs(total));
            local_min_gap = std::min(local_min_gap, gap);
            local_max_gap = std::max(local_max_gap, gap);
            local_min_magnitude = std::min(local_min_magnitude, magnitude);
            local_max_magnitude = std::max(local_max_magnitude, magnitude);
            local_min_eta_distance = std::min(local_min_eta_distance, std::abs(magnitude - lca_eta));
        }
        BranchMargins margins;
        margins.min_abs_total_density = pool_min(local_min_abs_density);
        margins.min_signed_saturation_gap = pool_min(local_min_gap);
        margins.max_signed_saturation_gap = pool_max(local_max_gap);
        margins.min_magnitude = pool_min(local_min_magnitude);
        margins.max_magnitude = pool_max(local_max_magnitude);
        margins.min_eta_distance = pool_min(local_min_eta_distance);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_BRANCH mode=" << mode
                      << " min_abs_total_density=" << margins.min_abs_total_density
                      << " min_signed_saturation_gap=" << margins.min_signed_saturation_gap
                      << " max_signed_saturation_gap=" << margins.max_signed_saturation_gap
                      << " min_magnitude=" << margins.min_magnitude << " max_magnitude=" << margins.max_magnitude
                      << " min_eta_distance=" << margins.min_eta_distance << '\n';
        }
        return margins;
    }

    void expect_vtxc_matches_returned_potential(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult result = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(result);
        ASSERT_EQ(potential.nr, 4);
        ASSERT_EQ(potential.nc, pw.nrxx);
        double local_inner_product = 0.0;
        for (int channel = 0; channel < 4; ++channel)
        {
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                local_inner_product += potential(channel, ir) * density[channel][ir];
            }
        }
        const double direct = pw.omega / pw.nxyz * pool_sum(local_inner_product);
        const double reported = std::get<1>(result);
        const double scale = std::max(1.0, std::max(std::abs(reported), std::abs(direct)));
        EXPECT_NEAR(reported, direct, 2.0e-12 * scale);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_VTXC mode=" << mode << " energy=" << std::get<0>(result)
                      << " reported=" << reported << " direct=" << direct
                      << " absolute_error=" << std::abs(reported - direct)
                      << " potential_fnv1a64=" << potential_fnv1a64(potential) << '\n';
        }
    }

    void expect_directional_derivatives_at_steps(const std::string& mode,
                                                 const Evaluator& evaluate,
                                                 const std::vector<double>& eps_values)
    {
        ASSERT_GE(eps_values.size(), 3U);
        const VxcResult base = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(base);
        ASSERT_EQ(potential.nr, 4);
        ASSERT_EQ(potential.nc, pw.nrxx);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_ENERGY mode=" << mode << " energy=" << std::get<0>(base)
                      << " vtxc=" << std::get<1>(base) << '\n';
        }

        for (int channel = 0; channel < 4; ++channel)
        {
            SCOPED_TRACE(std::string("density channel ") + std::to_string(channel));
            double local_analytic = 0.0;
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                local_analytic += potential(channel, ir) * perturbation[channel][ir];
            }
            const double analytic = pw.omega / pw.nxyz * pool_sum(local_analytic);
            std::vector<double> finite_difference(eps_values.size(), 0.0);
            std::vector<double> errors(eps_values.size(), 0.0);

            for (std::size_t ieps = 0; ieps < eps_values.size(); ++ieps)
            {
                const double eps = eps_values[ieps];
                for (int ir = 0; ir < pw.nrxx; ++ir)
                {
                    density[channel][ir] += eps * perturbation[channel][ir];
                }
                const double energy_plus = std::get<0>(evaluate());
                for (int ir = 0; ir < pw.nrxx; ++ir)
                {
                    density[channel][ir] -= 2.0 * eps * perturbation[channel][ir];
                }
                const double energy_minus = std::get<0>(evaluate());
                for (int ir = 0; ir < pw.nrxx; ++ir)
                {
                    density[channel][ir] += eps * perturbation[channel][ir];
                }
                finite_difference[ieps] = (energy_plus - energy_minus) / (2.0 * eps);
                errors[ieps] = std::abs(finite_difference[ieps] - analytic);
            }

            const double scale = std::max(1.0, std::max(std::abs(analytic), std::abs(finite_difference.back())));
            EXPECT_LE(*std::min_element(errors.begin(), errors.end()), 3.0e-8 * scale);
            EXPECT_LE(errors[1], 0.4 * errors[0] + 5.0e-9 * scale);
            EXPECT_LE(errors[2], 0.4 * errors[1] + 5.0e-9 * scale);
            if (is_pool_root())
            {
                for (std::size_t ieps = 0; ieps < eps_values.size(); ++ieps)
                {
                    const double relative_scale
                        = std::max(1.0e-30, std::max(std::abs(analytic), std::abs(finite_difference[ieps])));
                    const double order = ieps == 0 || errors[ieps] == 0.0
                                             ? 0.0
                                             : std::log(errors[ieps - 1] / errors[ieps]) / std::log(2.0);
                    std::cout << std::setprecision(17) << "NCGGA_FD mode=" << mode << " channel=" << channel
                              << " eps=" << eps_values[ieps] << " analytic=" << analytic
                              << " finite_difference=" << finite_difference[ieps] << " absolute_error=" << errors[ieps]
                              << " relative_error=" << errors[ieps] / relative_scale << " convergence_order=" << order
                              << '\n';
                }
            }
        }
    }

    void expect_core_directional_derivative(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult base = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(base);
        double local_analytic = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            local_analytic += potential(0, ir) * perturbation[0][ir];
        }
        const double analytic = pw.omega / pw.nxyz * pool_sum(local_analytic);
        const std::array<double, 4> eps_values = {{1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3}};
        std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};

        for (int ieps = 0; ieps < 4; ++ieps)
        {
            const double eps = eps_values[ieps];
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] += eps * perturbation[0][ir];
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double energy_plus = std::get<0>(evaluate());
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] -= 2.0 * eps * perturbation[0][ir];
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double energy_minus = std::get<0>(evaluate());
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] += eps * perturbation[0][ir];
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double finite_difference = (energy_plus - energy_minus) / (2.0 * eps);
            errors[ieps] = std::abs(finite_difference - analytic);
            if (is_pool_root())
            {
                std::cout << std::setprecision(17) << "NCGGA_CORE_FD mode=" << mode << " eps=" << eps
                          << " analytic=" << analytic << " finite_difference=" << finite_difference
                          << " absolute_error=" << errors[ieps] << '\n';
            }
        }
        const double scale = std::max(1.0, std::abs(analytic));
        EXPECT_LE(*std::min_element(errors.begin(), errors.end()), 3.0e-8 * scale);
        EXPECT_LE(errors[1], 0.4 * errors[0] + 5.0e-9 * scale);
        EXPECT_LE(errors[2], 0.4 * errors[1] + 5.0e-9 * scale);
    }

    void expect_core_translation_force(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult base = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(base);
        std::vector<std::complex<double>> reciprocal(pw.npw);
        std::vector<ModuleBase::Vector3<double>> core_gradient(pw.nrxx);
        pw.real2recip(core_density.data(), reciprocal.data());
        XC_Functional::grad_rho(reciprocal.data(), core_gradient.data(), &pw, pw.tpiba);
        double local_analytic = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            local_analytic += potential(0, ir) * core_gradient[ir].x;
        }
        const double analytic_force = pw.omega / pw.nxyz * pool_sum(local_analytic);
        const std::vector<double> original_core = core_density;
        const std::array<double, 4> eps_values = {{2.0e-2, 1.0e-2, 5.0e-3, 2.5e-3}};
        std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};
        for (int ieps = 0; ieps < 4; ++ieps)
        {
            const double eps = eps_values[ieps];
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] = original_core[ir] - eps * core_gradient[ir].x;
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double energy_plus = std::get<0>(evaluate());
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                core_density[ir] = original_core[ir] + eps * core_gradient[ir].x;
            }
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double energy_minus = std::get<0>(evaluate());
            std::copy(original_core.begin(), original_core.end(), core_density.begin());
            pw.real2recip(core_density.data(), core_density_reciprocal.data());
            const double finite_difference_force = -(energy_plus - energy_minus) / (2.0 * eps);
            errors[ieps] = std::abs(finite_difference_force - analytic_force);
            if (is_pool_root())
            {
                std::cout << std::setprecision(17) << "NCGGA_CORE_FORCE_FD mode=" << mode << " eps=" << eps
                          << " analytic=" << analytic_force << " finite_difference=" << finite_difference_force
                          << " absolute_error=" << errors[ieps] << '\n';
            }
        }
        const double scale = std::max(1.0, std::abs(analytic_force));
        EXPECT_LE(*std::min_element(errors.begin(), errors.end()), 3.0e-8 * scale);
        EXPECT_LE(errors[1], 0.4 * errors[0] + 5.0e-9 * scale);
        EXPECT_LE(errors[2], 0.4 * errors[1] + 5.0e-9 * scale);
    }

    void expect_core_repartition_invariance(const Evaluator& evaluate)
    {
        const VxcResult original = evaluate();
        const ModuleBase::matrix original_potential = std::get<2>(original);
        double local_expected_vtxc_change = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const double transfer = 0.04 * perturbation[0][ir];
            density[0][ir] += transfer;
            core_density[ir] -= transfer;
            local_expected_vtxc_change += original_potential(0, ir) * transfer;
        }
        pw.real2recip(core_density.data(), core_density_reciprocal.data());
        const VxcResult repartitioned = evaluate();
        const ModuleBase::matrix& repartitioned_potential = std::get<2>(repartitioned);
        EXPECT_NEAR(std::get<0>(original),
                    std::get<0>(repartitioned),
                    5.0e-11 * std::max(1.0, std::abs(std::get<0>(original))));
        for (int channel = 0; channel < 4; ++channel)
        {
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                EXPECT_NEAR(repartitioned_potential(channel, ir),
                            original_potential(channel, ir),
                            8.0e-11 * std::max(1.0, std::abs(original_potential(channel, ir))));
            }
        }
        const double expected_vtxc_change = pw.omega / pw.nxyz * pool_sum(local_expected_vtxc_change);
        const double actual_vtxc_change = std::get<1>(repartitioned) - std::get<1>(original);
        EXPECT_NEAR(actual_vtxc_change, expected_vtxc_change, 8.0e-11 * std::max(1.0, std::abs(expected_vtxc_change)));
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_CORE_REPARTITION expected_vtxc_change=" << expected_vtxc_change
                      << " actual_vtxc_change=" << actual_vtxc_change << '\n';
        }
    }

    void expect_local_rotation_torque(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult base = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(base);
        const std::vector<double> original_mx = density[1];
        const std::vector<double> original_my = density[2];
        double local_analytic = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            local_analytic
                += perturbation[0][ir] * (potential(2, ir) * original_mx[ir] - potential(1, ir) * original_my[ir]);
        }
        const double analytic = pw.omega / pw.nxyz * pool_sum(local_analytic);
        const std::array<double, 4> eps_values = {{1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3}};
        std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};

        for (int ieps = 0; ieps < 4; ++ieps)
        {
            const double eps = eps_values[ieps];
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                const double angle = eps * perturbation[0][ir];
                density[1][ir] = std::cos(angle) * original_mx[ir] - std::sin(angle) * original_my[ir];
                density[2][ir] = std::sin(angle) * original_mx[ir] + std::cos(angle) * original_my[ir];
            }
            const double energy_plus = std::get<0>(evaluate());
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                const double angle = -eps * perturbation[0][ir];
                density[1][ir] = std::cos(angle) * original_mx[ir] - std::sin(angle) * original_my[ir];
                density[2][ir] = std::sin(angle) * original_mx[ir] + std::cos(angle) * original_my[ir];
            }
            const double energy_minus = std::get<0>(evaluate());
            std::copy(original_mx.begin(), original_mx.end(), density[1].begin());
            std::copy(original_my.begin(), original_my.end(), density[2].begin());
            const double finite_difference = (energy_plus - energy_minus) / (2.0 * eps);
            errors[ieps] = std::abs(finite_difference - analytic);
            if (is_pool_root())
            {
                std::cout << std::setprecision(17) << "NCGGA_TORQUE_FD mode=" << mode << " eps=" << eps
                          << " analytic=" << analytic << " finite_difference=" << finite_difference
                          << " absolute_error=" << errors[ieps] << '\n';
            }
        }
        const double scale = std::max(1.0, std::abs(analytic));
        EXPECT_LE(*std::min_element(errors.begin(), errors.end()), 3.0e-8 * scale);
        EXPECT_LE(errors[1], 0.4 * errors[0] + 5.0e-9 * scale);
        EXPECT_LE(errors[2], 0.4 * errors[1] + 5.0e-9 * scale);
    }

    void expect_magnetization_inversion(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult original = evaluate();
        const ModuleBase::matrix original_potential = std::get<2>(original);
        for (int mu = 1; mu < 4; ++mu)
        {
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                density[mu][ir] = -density[mu][ir];
            }
        }
        const VxcResult inverted = evaluate();
        const ModuleBase::matrix& inverted_potential = std::get<2>(inverted);

        double local_charge_error = 0.0;
        double local_spin_error = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            local_charge_error
                = std::max(local_charge_error, std::abs(inverted_potential(0, ir) - original_potential(0, ir)));
            for (int mu = 1; mu < 4; ++mu)
            {
                local_spin_error
                    = std::max(local_spin_error, std::abs(inverted_potential(mu, ir) + original_potential(mu, ir)));
            }
        }
        const double charge_error = pool_max(local_charge_error);
        const double spin_error = pool_max(local_spin_error);
        const double energy_error = std::abs(std::get<0>(inverted) - std::get<0>(original));
        const double vtxc_error = std::abs(std::get<1>(inverted) - std::get<1>(original));
        const double energy_scale = std::max(1.0, std::abs(std::get<0>(original)));
        const double vtxc_scale = std::max(1.0, std::abs(std::get<1>(original)));
        EXPECT_LE(energy_error, 3.0e-11 * energy_scale);
        EXPECT_LE(vtxc_error, 3.0e-11 * vtxc_scale);
        EXPECT_LE(charge_error, 5.0e-11);
        EXPECT_LE(spin_error, 5.0e-11);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_INVERSION mode=" << mode << " energy_error=" << energy_error
                      << " vtxc_error=" << vtxc_error << " max_charge_error=" << charge_error
                      << " max_spin_error=" << spin_error << '\n';
        }
    }

    void expect_global_spin_rotation_covariance(const std::string& mode, const Evaluator& evaluate)
    {
        const VxcResult original = evaluate();
        const ModuleBase::matrix original_potential = std::get<2>(original);
        const double angle = 0.371;
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const double mx = density[1][ir];
            const double my = density[2][ir];
            density[1][ir] = cosine * mx - sine * my;
            density[2][ir] = sine * mx + cosine * my;
        }
        const VxcResult rotated = evaluate();
        const ModuleBase::matrix& rotated_potential = std::get<2>(rotated);

        double local_charge_error = 0.0;
        double local_spin_error = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            local_charge_error
                = std::max(local_charge_error, std::abs(rotated_potential(0, ir) - original_potential(0, ir)));
            const double expected_x = cosine * original_potential(1, ir) - sine * original_potential(2, ir);
            const double expected_y = sine * original_potential(1, ir) + cosine * original_potential(2, ir);
            local_spin_error = std::max(local_spin_error, std::abs(rotated_potential(1, ir) - expected_x));
            local_spin_error = std::max(local_spin_error, std::abs(rotated_potential(2, ir) - expected_y));
            local_spin_error
                = std::max(local_spin_error, std::abs(rotated_potential(3, ir) - original_potential(3, ir)));
        }
        const double charge_error = pool_max(local_charge_error);
        const double spin_error = pool_max(local_spin_error);
        const double energy_error = std::abs(std::get<0>(rotated) - std::get<0>(original));
        const double vtxc_error = std::abs(std::get<1>(rotated) - std::get<1>(original));
        const double energy_scale = std::max(1.0, std::abs(std::get<0>(original)));
        const double vtxc_scale = std::max(1.0, std::abs(std::get<1>(original)));
        EXPECT_LE(energy_error, 3.0e-11 * energy_scale);
        EXPECT_LE(vtxc_error, 3.0e-11 * vtxc_scale);
        EXPECT_LE(charge_error, 5.0e-11);
        EXPECT_LE(spin_error, 8.0e-11);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_GLOBAL_ROTATION mode=" << mode << " angle=" << angle
                      << " energy_error=" << energy_error << " vtxc_error=" << vtxc_error
                      << " max_charge_error=" << charge_error << " max_spin_error=" << spin_error << '\n';
        }
    }

    void expect_zero_magnetization_is_regular(const std::string& mode, const Evaluator& evaluate)
    {
        for (int mu = 1; mu < 4; ++mu)
        {
            std::fill(density[mu].begin(), density[mu].end(), 0.0);
        }
        const VxcResult result = evaluate();
        const ModuleBase::matrix& potential = std::get<2>(result);
        double local_maximum_spin_potential = 0.0;
        EXPECT_TRUE(std::isfinite(std::get<0>(result)));
        EXPECT_TRUE(std::isfinite(std::get<1>(result)));
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            EXPECT_TRUE(std::isfinite(potential(0, ir)));
            for (int mu = 1; mu < 4; ++mu)
            {
                EXPECT_TRUE(std::isfinite(potential(mu, ir)));
                local_maximum_spin_potential = std::max(local_maximum_spin_potential, std::abs(potential(mu, ir)));
            }
        }
        const double maximum_spin_potential = pool_max(local_maximum_spin_potential);
        EXPECT_DOUBLE_EQ(maximum_spin_potential, 0.0);
        if (is_pool_root())
        {
            std::cout << std::setprecision(17) << "NCGGA_ZERO_MAG mode=" << mode << " energy=" << std::get<0>(result)
                      << " vtxc=" << std::get<1>(result) << " max_spin_potential=" << maximum_spin_potential << '\n';
        }
    }
};

TEST_F(RealPwNcgga, PerturbationsSurviveThePwCutoff)
{
    for (int channel = 0; channel < 4; ++channel)
    {
        std::vector<std::complex<double>> reciprocal(pw.npw);
        pw.real2recip(perturbation[channel].data(), reciprocal.data());
        double local_norm2 = 0.0;
        for (int ig = 0; ig < pw.npw; ++ig)
        {
            local_norm2 += std::norm(reciprocal[ig]);
        }
        EXPECT_GT(pool_sum(local_norm2), 1.0e-4);
    }
}

TEST_F(RealPwNcgga, GradAndDivAreNegativeAdjoints)
{
    std::vector<double> scalar(pw.nrxx);
    std::vector<ModuleBase::Vector3<double>> vector_field(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        const int ix = ir / (pw.ny * pw.nplane);
        const int iy = (ir / pw.nplane) % pw.ny;
        const int iz = ir % pw.nplane + pw.startz_current;
        const double x = ModuleBase::TWO_PI * static_cast<double>(ix) / pw.nx;
        const double y = ModuleBase::TWO_PI * static_cast<double>(iy) / pw.ny;
        const double z = ModuleBase::TWO_PI * static_cast<double>(iz) / pw.nz;
        scalar[ir] = 0.4 * std::sin(3.0 * x + 0.2) - 0.3 * std::cos(8.0 * x - 0.1) + 0.2 * std::sin(y - z + 0.4);
        vector_field[ir].x = 0.7 * std::cos(2.0 * x + 0.3) + 0.2 * std::sin(7.0 * x);
        vector_field[ir].y = 0.3 * std::sin(y - 0.4) + 0.1 * std::cos(x + z);
        vector_field[ir].z = -0.25 * std::cos(z + 0.1) + 0.08 * std::sin(x - y);
    }

    std::vector<std::complex<double>> reciprocal(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<double> divergence(pw.nrxx);
    pw.real2recip(scalar.data(), reciprocal.data());
    XC_Functional::grad_rho(reciprocal.data(), gradient.data(), &pw, pw.tpiba);
    XC_Functional::grad_dot(vector_field.data(), divergence.data(), &pw, pw.tpiba);

    double local_identity = 0.0;
    double local_norm = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        const double left = vector_field[ir] * gradient[ir];
        const double right = divergence[ir] * scalar[ir];
        local_identity += left + right;
        local_norm += std::abs(left) + std::abs(right);
    }
    const double identity = pw.omega / pw.nxyz * pool_sum(local_identity);
    const double norm = pw.omega / pw.nxyz * pool_sum(local_norm);
    EXPECT_GT(norm, 1.0e-4);
    EXPECT_LE(std::abs(identity), 5.0e-11 * std::max(1.0, norm));
    if (is_pool_root())
    {
        std::cout << std::setprecision(17) << "NCGGA_ADJOINT identity=" << identity << " norm=" << norm
                  << " scaled_error=" << std::abs(identity) / norm << '\n';
    }
}

TEST_F(RealPwNcgga, ProjectedLcaGraphRequiresDiscreteFluxReverse)
{
    // This eta is the documented gga_grad=2 policy.  Keep it local so the
    // identical test-only patch compiles on the behavior commit's parent.
    constexpr double lca_eta = 1.0e-3;
    typedef std::array<std::vector<ModuleBase::Vector3<double>>, 3> Gradients;

    const auto gradients = [&]() {
        Gradients result;
        std::vector<std::complex<double>> reciprocal(pw.npw);
        for (int mu = 0; mu < 3; ++mu)
        {
            result[mu].resize(pw.nrxx);
            pw.real2recip(density[mu + 1].data(), reciprocal.data());
            XC_Functional::grad_rho(reciprocal.data(), result[mu].data(), &pw, pw.tpiba);
        }
        return result;
    };

    const auto graph_energy = [&]() {
        const Gradients grad_m = gradients();
        double local_energy = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const std::array<double, 3> magnetization = {{density[1][ir], density[2][ir], density[3][ir]}};
            const ModuleXC::NcggaRadialPoint radial = ModuleXC::make_ncgga_radial_point(magnetization, lca_eta);
            ModuleBase::Vector3<double> projected;
            for (int nu = 0; nu < 3; ++nu)
            {
                projected += radial.gradient[nu] * grad_m[nu][ir];
            }
            local_energy += 0.5 * (projected * projected);
        }
        return pw.omega / pw.nxyz * pool_sum(local_energy);
    };

    const Gradients grad_m = gradients();
    std::vector<ModuleBase::Vector3<double>> projected(pw.nrxx);
    std::array<std::vector<double>, 3> exact_potential;
    std::array<std::vector<double>, 3> old_surrogate;
    std::vector<double> projected_divergence(pw.nrxx);
    std::vector<ModuleBase::Vector3<double>> flux(pw.nrxx);
    std::vector<double> divergence(pw.nrxx);
    for (int mu = 0; mu < 3; ++mu)
    {
        exact_potential[mu].resize(pw.nrxx);
        old_surrogate[mu].resize(pw.nrxx);
    }
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        const std::array<double, 3> magnetization = {{density[1][ir], density[2][ir], density[3][ir]}};
        const ModuleXC::NcggaRadialPoint radial = ModuleXC::make_ncgga_radial_point(magnetization, lca_eta);
        for (int nu = 0; nu < 3; ++nu)
        {
            projected[ir] += radial.gradient[nu] * grad_m[nu][ir];
        }
    }
    XC_Functional::grad_dot(projected.data(), projected_divergence.data(), &pw, pw.tpiba);

    double maximum_surrogate_error = 0.0;
    double maximum_exact_error = 0.0;
    const std::array<double, 5> eps_values = {{1.0e-3, 5.0e-4, 2.5e-4, 1.25e-4, 6.25e-5}};
    for (int mu = 0; mu < 3; ++mu)
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const std::array<double, 3> magnetization = {{density[1][ir], density[2][ir], density[3][ir]}};
            const ModuleXC::NcggaRadialPoint radial = ModuleXC::make_ncgga_radial_point(magnetization, lca_eta);
            flux[ir] = radial.gradient[mu] * projected[ir];
        }
        XC_Functional::grad_dot(flux.data(), divergence.data(), &pw, pw.tpiba);

        double local_exact_projection = 0.0;
        double local_surrogate_projection = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const std::array<double, 3> magnetization = {{density[1][ir], density[2][ir], density[3][ir]}};
            const ModuleXC::NcggaRadialPoint radial = ModuleXC::make_ncgga_radial_point(magnetization, lca_eta);
            double local_response = 0.0;
            for (int nu = 0; nu < 3; ++nu)
            {
                local_response += radial.jacobian(nu, mu) * (projected[ir] * grad_m[nu][ir]);
            }
            exact_potential[mu][ir] = local_response - divergence[ir];
            old_surrogate[mu][ir] = -radial.gradient[mu] * projected_divergence[ir];
            local_exact_projection += exact_potential[mu][ir] * perturbation[mu + 1][ir];
            local_surrogate_projection += old_surrogate[mu][ir] * perturbation[mu + 1][ir];
        }
        const double exact = pw.omega / pw.nxyz * pool_sum(local_exact_projection);
        const double surrogate = pw.omega / pw.nxyz * pool_sum(local_surrogate_projection);

        std::array<double, 5> exact_errors = {{0.0, 0.0, 0.0, 0.0, 0.0}};
        std::array<double, 5> surrogate_errors = {{0.0, 0.0, 0.0, 0.0, 0.0}};
        for (int ieps = 0; ieps < 5; ++ieps)
        {
            const double eps = eps_values[ieps];
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                density[mu + 1][ir] += eps * perturbation[mu + 1][ir];
            }
            const double energy_plus = graph_energy();
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                density[mu + 1][ir] -= 2.0 * eps * perturbation[mu + 1][ir];
            }
            const double energy_minus = graph_energy();
            for (int ir = 0; ir < pw.nrxx; ++ir)
            {
                density[mu + 1][ir] += eps * perturbation[mu + 1][ir];
            }
            const double finite_difference = (energy_plus - energy_minus) / (2.0 * eps);
            exact_errors[ieps] = std::abs(exact - finite_difference);
            surrogate_errors[ieps] = std::abs(surrogate - finite_difference);
            maximum_exact_error = std::max(maximum_exact_error, exact_errors[ieps]);
            maximum_surrogate_error = std::max(maximum_surrogate_error, surrogate_errors[ieps]);
            if (is_pool_root())
            {
                std::cout << std::setprecision(17) << "NCGGA_PROJECTED_REVERSE channel=" << mu + 1 << " eps=" << eps
                          << " exact=" << exact << " old_surrogate=" << surrogate
                          << " finite_difference=" << finite_difference << " exact_error=" << exact_errors[ieps]
                          << " surrogate_error=" << surrogate_errors[ieps] << '\n';
            }
        }
        const double scale = std::max(1.0, std::max(std::abs(exact), std::abs(surrogate)));
        EXPECT_LE(*std::min_element(exact_errors.begin(), exact_errors.end()), 2.0e-8 * scale);
        EXPECT_GT(*std::min_element(surrogate_errors.begin(), surrogate_errors.end()), 2.0e-7 * scale);
    }
    EXPECT_GT(maximum_surrogate_error, 100.0 * maximum_exact_error);
}

#ifdef __LIBXC

TEST_F(RealPwNcgga, LibxcGgaGrad2VtxcBookkeepingUsesFinalReturnedPotential)
{
    const std::vector<int> functionals = {XC_LDA_X, XC_GGA_C_PBE};
    const Evaluator evaluate = [this, &functionals]() { return evaluate_libxc(functionals, nullptr); };
    expect_vtxc_matches_returned_potential("libxc_gga2_mixed", evaluate);
}

TEST_F(RealPwNcgga, LibxcGgaGrad2VtxcEqualsFinalValencePotentialInnerProduct)
{
    report_branch_margins("libxc_vtxc_smooth");
    expect_vtxc_matches_returned_potential("libxc_gga2_smooth", [this]() { return evaluate_libxc_gga(); });

    set_negative_gga_state();
    report_branch_margins("libxc_vtxc_negative");
    expect_vtxc_matches_returned_potential("libxc_gga2_negative", [this]() { return evaluate_libxc_gga(); });

    set_saturated_gga_state();
    report_branch_margins("libxc_vtxc_saturated");
    expect_vtxc_matches_returned_potential("libxc_gga2_saturated", [this]() { return evaluate_libxc_gga(); });

    set_inside_eta_state();
    report_branch_margins("libxc_vtxc_inside_eta");
    expect_vtxc_matches_returned_potential("libxc_gga2_inside_eta", [this]() { return evaluate_libxc_gga(); });
}

TEST_F(RealPwNcgga, LibxcGgaGrad2IsDiscreteGradientOnSmoothProjectedBranch)
{
    const BranchMargins margins = report_branch_margins("libxc_smooth");
    EXPECT_GT(margins.min_abs_total_density, 1.0);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    EXPECT_GT(margins.min_eta_distance, 0.3);
    expect_directional_derivatives_at_steps("libxc_gga2_smooth",
                                            [this]() { return evaluate_libxc_gga(); },
                                            {1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, LibxcGgaGrad2IsDiscreteGradientInsideRadialEta)
{
    set_inside_eta_state();
    const BranchMargins margins = report_branch_margins("libxc_inside_eta");
    EXPECT_GT(margins.min_abs_total_density, 0.025);
    EXPECT_LT(margins.max_magnitude, 6.0e-4);
    EXPECT_GT(margins.min_eta_distance, 4.0e-4);
    expect_directional_derivatives_at_steps(
        "libxc_gga2_inside_eta",
        [this]() { return evaluate_libxc_gga(); },
        {2.0e-4, 1.0e-4, 5.0e-5, 2.5e-5, 1.25e-5, 6.25e-6, 3.125e-6, 1.5625e-6, 7.8125e-7, 3.90625e-7});
}

TEST_F(RealPwNcgga, LibxcGgaGrad2DifferentiatesNegativeDensityBranch)
{
    set_negative_gga_state();
    const BranchMargins margins = report_branch_margins("libxc_negative");
    EXPECT_GT(margins.min_abs_total_density, 1.3);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    expect_directional_derivatives_at_steps("libxc_gga2_negative",
                                            [this]() { return evaluate_libxc_gga(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, LibxcGgaGrad2DifferentiatesSaturatedGgaBranch)
{
    set_saturated_gga_state();
    const BranchMargins margins = report_branch_margins("libxc_saturated");
    EXPECT_LT(margins.max_signed_saturation_gap, -0.1);
    expect_directional_derivatives_at_steps("libxc_gga2_saturated",
                                            [this]() { return evaluate_libxc_gga(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, LibxcLdaGgaGrad2DifferentiatesLocalMapBranches)
{
    set_uniform_state(-1.4, {{0.20, -0.16, 0.18}});
    const BranchMargins negative = report_branch_margins("libxc_lda_negative");
    EXPECT_GT(negative.min_abs_total_density, 1.3);
    EXPECT_GT(negative.min_signed_saturation_gap, 1.0);
    expect_directional_derivatives_at_steps("libxc_lda_gga2_negative",
                                            [this]() { return evaluate_libxc_lda(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});

    set_uniform_state(0.45, {{0.65, 0.30, 0.20}});
    const BranchMargins saturated = report_branch_margins("libxc_lda_saturated");
    EXPECT_LT(saturated.max_signed_saturation_gap, -0.25);
    expect_directional_derivatives_at_steps("libxc_lda_gga2_saturated",
                                            [this]() { return evaluate_libxc_lda(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, LibxcGgaGrad2DifferentiatesCoreDensityAndLocalRotation)
{
    const BranchMargins margins = report_branch_margins("libxc_core_rotation");
    EXPECT_GT(margins.min_abs_total_density, 1.0);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    expect_core_directional_derivative("libxc_gga2_core", [this]() { return evaluate_libxc_gga(); });
    expect_core_translation_force("libxc_gga2_core_translation", [this]() { return evaluate_libxc_gga(); });
    expect_local_rotation_torque("libxc_gga2_local_rotation", [this]() { return evaluate_libxc_gga(); });
    expect_core_repartition_invariance([this]() { return evaluate_libxc_gga(); });
}

TEST_F(RealPwNcgga, LibxcGgaGrad2RespectsSpinInversionAndZeroLimit)
{
    report_branch_margins("libxc_symmetry");
    expect_magnetization_inversion("libxc_gga2", [this]() { return evaluate_libxc_gga(); });
    expect_zero_magnetization_is_regular("libxc_gga2", [this]() { return evaluate_libxc_gga(); });
}

TEST_F(RealPwNcgga, LibxcGgaGrad2IsCovariantUnderGlobalSpinRotation)
{
    report_branch_margins("libxc_global_rotation");
    expect_global_spin_rotation_covariance("libxc_gga2", [this]() { return evaluate_libxc_gga(); });
}

TEST_F(RealPwNcgga, LibxcGgaGrad2DifferentiatesNonuniformFunctionalScaling)
{
    const std::map<int, double> scaling_factor = {{XC_GGA_X_PBE, 0.37}, {XC_GGA_C_PBE, 1.23}};
    const BranchMargins margins = report_branch_margins("libxc_scaled");
    EXPECT_GT(margins.min_abs_total_density, 1.0);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    if (is_pool_root())
    {
        std::cout << std::setprecision(17) << "LIBXC_SCALING mode=libxc_gga2_scaled"
                  << " exchange=" << scaling_factor.at(XC_GGA_X_PBE)
                  << " correlation=" << scaling_factor.at(XC_GGA_C_PBE) << '\n';
    }
    const Evaluator evaluate = [this, &scaling_factor]() { return evaluate_libxc_gga(&scaling_factor); };

    // A finite difference alone cannot prove that component scaling was
    // applied: energy and potential could omit the same factors and remain
    // self-consistent.  Independently require linearity against unscaled
    // exchange-only and correlation-only evaluations.
    const VxcResult exchange = evaluate_libxc({XC_GGA_X_PBE}, nullptr);
    const VxcResult correlation = evaluate_libxc({XC_GGA_C_PBE}, nullptr);
    const VxcResult scaled = evaluate();
    const double exchange_factor = scaling_factor.at(XC_GGA_X_PBE);
    const double correlation_factor = scaling_factor.at(XC_GGA_C_PBE);
    const double expected_energy
        = exchange_factor * std::get<0>(exchange) + correlation_factor * std::get<0>(correlation);
    const double expected_vtxc
        = exchange_factor * std::get<1>(exchange) + correlation_factor * std::get<1>(correlation);
    const double energy_error = std::abs(std::get<0>(scaled) - expected_energy);
    const double vtxc_error = std::abs(std::get<1>(scaled) - expected_vtxc);
    double local_potential_error = 0.0;
    for (int channel = 0; channel < 4; ++channel)
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const double expected = exchange_factor * std::get<2>(exchange)(channel, ir)
                                    + correlation_factor * std::get<2>(correlation)(channel, ir);
            local_potential_error
                = std::max(local_potential_error, std::abs(std::get<2>(scaled)(channel, ir) - expected));
        }
    }
    const double potential_error = pool_max(local_potential_error);
    EXPECT_LE(energy_error, 5.0e-11 * std::max(1.0, std::abs(expected_energy)));
    EXPECT_LE(vtxc_error, 5.0e-11 * std::max(1.0, std::abs(expected_vtxc)));
    EXPECT_LE(potential_error, 8.0e-11);
    if (is_pool_root())
    {
        std::cout << std::setprecision(17) << "LIBXC_SCALING_LINEARITY mode=libxc_gga2_scaled"
                  << " energy_error=" << energy_error << " vtxc_error=" << vtxc_error
                  << " max_potential_error=" << potential_error << '\n';
    }

    expect_directional_derivatives_at_steps("libxc_gga2_scaled", evaluate, {1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
    expect_vtxc_matches_returned_potential("libxc_gga2_scaled", evaluate);
}

TEST_F(RealPwNcgga, LibxcGgaGrad2DifferentiatesMixedLdaAndGgaComponents)
{
    const std::vector<int> functionals = {XC_LDA_X, XC_GGA_C_PBE};
    const Evaluator evaluate = [this, &functionals]() { return evaluate_libxc(functionals, nullptr); };
    const BranchMargins margins = report_branch_margins("libxc_mixed_lda_gga");
    EXPECT_GT(margins.min_abs_total_density, 1.0);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    expect_directional_derivatives_at_steps("libxc_gga2_mixed_lda_gga",
                                            evaluate,
                                            {1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
    expect_vtxc_matches_returned_potential("libxc_gga2_mixed_lda_gga", evaluate);
}








#endif

TEST_F(RealPwNcgga, BuiltinGgaGrad2VtxcEqualsFinalValencePotentialInnerProduct)
{
    expect_vtxc_matches_returned_potential("smooth", [this]() { return evaluate_builtin(); });

    set_negative_gga_state();
    expect_vtxc_matches_returned_potential("negative", [this]() { return evaluate_builtin(); });

    set_saturated_gga_state();
    expect_vtxc_matches_returned_potential("saturated", [this]() { return evaluate_builtin(); });

    set_inside_eta_state();
    expect_vtxc_matches_returned_potential("inside_eta", [this]() { return evaluate_builtin(); });
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2IsDiscreteGradientOnSmoothProjectedBranch)
{
    const BranchMargins margins = report_branch_margins("smooth");
    EXPECT_GT(margins.min_abs_total_density, 1.0);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    EXPECT_GT(margins.min_eta_distance, 0.3);
    expect_directional_derivatives_at_steps("builtin_gga2_smooth",
                                            [this]() { return evaluate_builtin(); },
                                            {1.0e-2, 5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2IsDiscreteGradientInsideRadialEta)
{
    set_inside_eta_state();
    const BranchMargins margins = report_branch_margins("inside_eta");
    EXPECT_GT(margins.min_abs_total_density, 0.025);
    EXPECT_LT(margins.max_magnitude, 6.0e-4);
    EXPECT_GT(margins.min_eta_distance, 4.0e-4);
    expect_directional_derivatives_at_steps(
        "builtin_gga2_inside_eta",
        [this]() { return evaluate_builtin(); },
        {2.0e-4, 1.0e-4, 5.0e-5, 2.5e-5, 1.25e-5, 6.25e-6, 3.125e-6, 1.5625e-6, 7.8125e-7, 3.90625e-7});
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2DifferentiatesNegativeDensityBranch)
{
    set_negative_gga_state();
    const BranchMargins margins = report_branch_margins("negative");
    EXPECT_GT(margins.min_abs_total_density, 1.3);
    EXPECT_GT(margins.min_signed_saturation_gap, 0.8);
    expect_directional_derivatives_at_steps("builtin_gga2_negative",
                                            [this]() { return evaluate_builtin(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2DifferentiatesSaturatedGgaBranch)
{
    set_saturated_gga_state();
    const BranchMargins margins = report_branch_margins("saturated");
    EXPECT_LT(margins.max_signed_saturation_gap, -0.1);
    expect_directional_derivatives_at_steps("builtin_gga2_saturated_gga",
                                            [this]() { return evaluate_builtin(); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, BuiltinLdaGgaGrad2DifferentiatesLocalMapBranches)
{
    set_uniform_state(-1.4, {{0.20, -0.16, 0.18}});
    expect_directional_derivatives_at_steps("builtin_lda_gga2_negative",
                                            [this]() { return evaluate_builtin("PZ"); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});

    set_uniform_state(0.45, {{0.65, 0.30, 0.20}});
    expect_directional_derivatives_at_steps("builtin_lda_gga2_saturated",
                                            [this]() { return evaluate_builtin("PZ"); },
                                            {5.0e-3, 2.5e-3, 1.25e-3, 6.25e-4});
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2DifferentiatesCoreDensityAndLocalRotation)
{
    expect_core_directional_derivative("builtin_gga2_core", [this]() { return evaluate_builtin(); });
    expect_core_translation_force("builtin_gga2_core_translation", [this]() { return evaluate_builtin(); });
    expect_local_rotation_torque("builtin_gga2_rotation", [this]() { return evaluate_builtin(); });
    expect_core_repartition_invariance([this]() { return evaluate_builtin(); });
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2RespectsSpinSymmetriesAndZeroLimit)
{
    const VxcResult original = evaluate_builtin();
    const ModuleBase::matrix original_potential = std::get<2>(original);
    for (int mu = 1; mu < 4; ++mu)
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            density[mu][ir] = -density[mu][ir];
        }
    }
    const VxcResult inverted = evaluate_builtin();
    const ModuleBase::matrix& inverted_potential = std::get<2>(inverted);
    EXPECT_NEAR(std::get<0>(original), std::get<0>(inverted), 3.0e-11 * std::max(1.0, std::abs(std::get<0>(original))));
    EXPECT_NEAR(std::get<1>(original), std::get<1>(inverted), 3.0e-11 * std::max(1.0, std::abs(std::get<1>(original))));
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        EXPECT_NEAR(inverted_potential(0, ir), original_potential(0, ir), 3.0e-11);
        for (int mu = 1; mu < 4; ++mu)
        {
            EXPECT_NEAR(inverted_potential(mu, ir),
                        -original_potential(mu, ir),
                        5.0e-11 * std::max(1.0, std::abs(original_potential(mu, ir))));
        }
    }

    for (int mu = 1; mu < 4; ++mu)
    {
        std::fill(density[mu].begin(), density[mu].end(), 0.0);
    }
    const VxcResult zero = evaluate_builtin();
    EXPECT_TRUE(std::isfinite(std::get<0>(zero)));
    EXPECT_TRUE(std::isfinite(std::get<1>(zero)));
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        EXPECT_TRUE(std::isfinite(std::get<2>(zero)(0, ir)));
        EXPECT_DOUBLE_EQ(std::get<2>(zero)(1, ir), 0.0);
        EXPECT_DOUBLE_EQ(std::get<2>(zero)(2, ir), 0.0);
        EXPECT_DOUBLE_EQ(std::get<2>(zero)(3, ir), 0.0);
    }
}

TEST_F(RealPwNcgga, BuiltinGgaGrad2IsCovariantUnderGlobalSpinRotation)
{
    const VxcResult original = evaluate_builtin();
    const ModuleBase::matrix original_potential = std::get<2>(original);
    const double angle = 0.371;
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        const double mx = density[1][ir];
        const double my = density[2][ir];
        density[1][ir] = cosine * mx - sine * my;
        density[2][ir] = sine * mx + cosine * my;
    }
    const VxcResult rotated = evaluate_builtin();
    const ModuleBase::matrix& rotated_potential = std::get<2>(rotated);
    EXPECT_NEAR(std::get<0>(original), std::get<0>(rotated), 3.0e-11 * std::max(1.0, std::abs(std::get<0>(original))));
    EXPECT_NEAR(std::get<1>(original), std::get<1>(rotated), 3.0e-11 * std::max(1.0, std::abs(std::get<1>(original))));
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        EXPECT_NEAR(rotated_potential(0, ir), original_potential(0, ir), 3.0e-11);
        const double expected_x = cosine * original_potential(1, ir) - sine * original_potential(2, ir);
        const double expected_y = sine * original_potential(1, ir) + cosine * original_potential(2, ir);
        EXPECT_NEAR(rotated_potential(1, ir), expected_x, 5.0e-11 * std::max(1.0, std::abs(expected_x)));
        EXPECT_NEAR(rotated_potential(2, ir), expected_y, 5.0e-11 * std::max(1.0, std::abs(expected_y)));
        EXPECT_NEAR(rotated_potential(3, ir),
                    original_potential(3, ir),
                    5.0e-11 * std::max(1.0, std::abs(original_potential(3, ir))));
    }
}









} // namespace

#ifdef __MPI
int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, test_size, threads, test_rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
#endif
