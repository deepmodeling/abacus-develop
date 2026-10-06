#include "sccs_response.h"
#include "sccs_parameters.h"
#include "sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <sstream>
#include <utility>

namespace ModuleSccs
{
namespace
{
void prepare_cavity(const std::vector<double>& density,
                    const CavityParameters& cavity,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    SccsResponse& response,
                    std::vector<double>& coefficient)
{
    const std::size_t size = density.size();
    response.solute.resize(size);
    response.dsolute_drho.resize(size);
    response.epsilon.resize(size);
    response.depsilon_drho.resize(size);
    response.grad_log_epsilon.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const CavityPoint point = evaluate_cavity(density[i], cavity);
        response.solute[i] = point.solute;
        response.dsolute_drho[i] = point.dsolute_drho;
        response.epsilon[i] = point.epsilon;
        response.depsilon_drho[i] = point.depsilon_drho;
    }

    std::vector<std::complex<double>> density_g(basis.npw);
    basis.real2recip(density.data(), density_g.data());
    std::vector<ModuleBase::Vector3<double>> gradient(size);
    std::vector<double> laplacian(size);
    // These XC helpers use G in units of tpiba and ABACUS's normalized FFT.
    XC_Functional::grad_rho(density_g.data(), gradient.data(), &basis, tpiba);
    XC_Functional::laplacian_rho(density_g.data(), laplacian.data(), &basis, tpiba);
    const double density_ratio = cavity.density_max / cavity.density_min;
    const double width = std::log(density_ratio);
    const double log_bulk = std::log(cavity.epsilon_bulk);
    coefficient.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        double second_log = 0.0;
        if (density[i] > cavity.density_min && density[i] < cavity.density_max)
        {
            const double local_ratio = cavity.density_max / density[i];
            const double x = std::log(local_ratio) / width;
            const double angle = ModuleBase::TWO_PI * x;
            second_log = log_bulk * (1.0 - std::cos(angle) + ModuleBase::TWO_PI * std::sin(angle) / width)
                         / (width * density[i] * density[i]);
        }
        const double first_log = response.depsilon_drho[i] / response.epsilon[i];
        const ModuleBase::Vector3<double>& local_gradient = gradient[i];
        const double gradient_square = local_gradient * local_gradient;
        response.grad_log_epsilon[i] = local_gradient * first_log;
        const double lap_log = first_log * laplacian[i] + second_log * gradient_square;
        coefficient[i] = response.epsilon[i] * (0.5 * lap_log + 0.25 * first_log * first_log * gradient_square)
                         / ModuleBase::FOUR_PI;
    }
}

// Pool-reduced norms, so every rank takes the same convergence decision. A
// non-finite residual fails convergence and stops at the next CG step.
void residual_norms(const std::vector<double>& values,
                    const ModulePW::PW_Basis& basis,
                    double& rms,
                    double& maximum)
{
    double square = 0.0;
    maximum = 0.0;
    for (double value : values)
    {
        square += value * value;
        const double magnitude = std::abs(value);
        maximum = std::max(maximum, magnitude);
    }
    Parallel_Reduce::reduce_pool(square);
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, maximum);
    const double mean_square = square / basis.nxyz;
    rms = std::sqrt(mean_square);
}

bool converged(const PolarizationResult& result, const PolarizationSolverParameters& solver)
{
    return result.residual_rms <= solver.tolerance_rms && result.residual_max <= solver.tolerance_max;
}

double grid_dot(const std::vector<double>& left,
                const std::vector<double>& right,
                const ModulePW::PW_Basis& basis)
{
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        value += left[i] * right[i];
    }
    Parallel_Reduce::reduce_pool(value);
    return value * basis.omega / basis.nxyz;
}

// P r = eps^-1/2 G eps^-1/2 r; only FFT scratch survives an application.
void apply_preconditioner(const std::vector<double>& rhs,
                          const std::vector<double>& invsqrt,
                          PeriodicCoulombOperator& coulomb,
                          std::vector<double>& weighted,
                          std::vector<double>& value)
{
    const std::size_t size = rhs.size();
    weighted.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        weighted[i] = rhs[i] * invsqrt[i];
    }
    coulomb.apply_potential(weighted, value);
    for (std::size_t i = 0; i < size; ++i)
    {
        value[i] *= invsqrt[i];
    }
}
} // namespace

void solve_sccs_response(const std::vector<double>& density,
                         const std::vector<double>& charge,
                         const CavityParameters& cavity,
                         const PolarizationSolverParameters& solver,
                         const std::vector<double>& initial_potential,
                         const ModulePW::PW_Basis& basis,
                         double tpiba,
                         SccsResponse& result)
{
    // Determine start mode collectively, including ranks with no real-space points.
    double initial_count = initial_potential.size();
    Parallel_Reduce::reduce_pool(initial_count);
    const bool warm_start = initial_count != 0.0;

    SccsResponse candidate;
    std::vector<double> coefficient;
    prepare_cavity(density, cavity, basis, tpiba, candidate, coefficient);
    const std::size_t size = density.size();
    std::vector<double> invsqrt(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        invsqrt[i] = 1.0 / std::sqrt(candidate.epsilon[i]);
    }
    PeriodicCoulombOperator coulomb(basis, tpiba);
    std::vector<double> residual = charge;
    std::vector<double> potential(size, 0.0);
    std::vector<double> direction(size, 0.0);
    std::vector<double> image(size, 0.0);
    std::vector<double> weighted;
    std::vector<double> z;
    PolarizationResult& polarization = candidate.polarization;
    residual_norms(residual, basis, polarization.residual_rms, polarization.residual_max);
    // Retain the old fixed-point warm start only when it reduces the charge residual.
    if (!converged(polarization, solver) && warm_start)
    {
        std::vector<double> guess_residual(size);
        for (std::size_t i = 0; i < size; ++i)
        {
            guess_residual[i] = charge[i] - coefficient[i] * initial_potential[i];
        }
        apply_preconditioner(guess_residual, invsqrt, coulomb, weighted, z);
        for (std::size_t i = 0; i < size; ++i)
        {
            guess_residual[i] = coefficient[i] * (initial_potential[i] - z[i]);
        }
        double guess_rms = 0.0;
        double guess_max = 0.0;
        residual_norms(guess_residual, basis, guess_rms, guess_max);
        if (guess_rms < polarization.residual_rms)
        {
            potential.swap(z);
            residual.swap(guess_residual);
            polarization.warm_started = true;
            polarization.residual_rms = guess_rms;
            polarization.residual_max = guess_max;
        }
    }
    double old_rz = 0.0;
    for (int iteration = 1; !converged(polarization, solver) && iteration <= solver.max_iterations; ++iteration)
    {
        apply_preconditioner(residual, invsqrt, coulomb, weighted, z);
        const double rz = grid_dot(residual, z, basis);
        if (!std::isfinite(rz) || std::abs(rz) < 1e-30)
        {
            ModuleBase::WARNING_QUIT("ModuleSccs::solve_sccs_response",
                                     "SCCS sqrt-CG has a null or nonfinite preconditioned residual");
        }
        const double beta = std::abs(old_rz) > 1e-30 ? rz / old_rz : 0.0;
        old_rz = rz;
        for (std::size_t i = 0; i < size; ++i)
        {
            direction[i] = z[i] + beta * direction[i];
            image[i] = coefficient[i] * z[i] + residual[i] + beta * image[i];
        }
        const double curvature = grid_dot(direction, image, basis);
        if (!std::isfinite(curvature) || curvature == 0.0)
        {
            ModuleBase::WARNING_QUIT("ModuleSccs::solve_sccs_response", "SCCS sqrt-CG has invalid curvature");
        }
        const double alpha = rz / curvature;
        for (std::size_t i = 0; i < size; ++i)
        {
            potential[i] += alpha * direction[i];
            residual[i] -= alpha * image[i];
        }
        polarization.iterations = iteration;
        residual_norms(residual, basis, polarization.residual_rms, polarization.residual_max);
    }
    if (!converged(polarization, solver))
    {
        std::ostringstream message;
        message << "SCCS sqrt-CG did not reach both residual tolerances within sccs_maxiter = "
                << solver.max_iterations << " iterations (RMS " << polarization.residual_rms << ", MAX "
                << polarization.residual_max << ")";
        ModuleBase::WARNING_QUIT("ModuleSccs::solve_sccs_response", message.str());
    }

    candidate.restart_potential = potential;
    double mean = 0.0;
    for (double value : potential)
    {
        mean += value;
    }
    Parallel_Reduce::reduce_pool(mean);
    mean /= basis.nxyz;
    for (double& value : potential)
    {
        value -= mean;
    }
    std::vector<std::complex<double>> potential_g(basis.npw);
    basis.real2recip(potential.data(), potential_g.data());
    polarization.gradient.resize(size);
    XC_Functional::grad_rho(potential_g.data(), polarization.gradient.data(), &basis, tpiba);
    polarization.potential.swap(potential);
    candidate.cavity_potential.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const ModuleBase::Vector3<double>& gradient = polarization.gradient[i];
        const double gradient_square = gradient * gradient;
        candidate.cavity_potential[i] = -candidate.depsilon_drho[i] * gradient_square / (8.0 * ModuleBase::PI);
    }
    result = std::move(candidate);
}
} // namespace ModuleSccs
