#include "sccs_response.h"
#include "sccs_cav_deriv.h"
#include "sccs_lowpass.h"
#include "sccs_parameters.h"
#include "sccs_pw_coulomb.h"
#include "sccs_thread_sum.h"

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
// sqrt(eps) v = C_PCC(s), s = (q - f v)/sqrt(eps). The far-field
// screening charge is int(s)/sqrt(eps_bulk) - int(q), including its physical sign.
void finish_open_bc(const std::vector<double>& charge,
                    const std::vector<double>& coefficient,
                    const std::vector<double>& invsqrt,
                    const CavityParameters& cavity,
                    const ModulePW::PW_Basis& basis,
                    SccsResponse& result)
{
    const std::size_t size = charge.size();
    const std::vector<double>& potential = result.polarization.potential;
    const auto add_point = [&](std::size_t i, std::array<double, 2>& sums) {
        sums[0] += (charge[i] - coefficient[i] * potential[i]) * invsqrt[i];
        sums[1] += charge[i];
    };
    const std::array<double, 2> sums = thread_sums<2>(size, add_point);
    double source_sum = sums[0];
    double solute_sum = sums[1];
    Parallel_Reduce::reduce_pool(source_sum);
    Parallel_Reduce::reduce_pool(solute_sum);
    const double volume_element = basis.omega / basis.nxyz;
    const double bulk_invsqrt = 1.0 / std::sqrt(cavity.epsilon_bulk);
    result.far_field_polarization_charge = (source_sum * bulk_invsqrt - solute_sum) * volume_element;
}

// Pool-reduced norms, so every rank takes the same convergence decision. A
// non-finite residual fails convergence and stops at the next CG step.
void residual_norms(const std::vector<double>& values,
                    const ModulePW::PW_Basis& basis,
                    double& rms,
                    double& maximum)
{
    const auto add_square = [&values](std::size_t i, std::array<double, 1>& sum) {
        sum[0] += values[i] * values[i];
    };
    double square = thread_sums<1>(values.size(), add_square)[0];
    double largest = 0.0;
    const long size = values.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max : largest)
#endif
    for (long i = 0; i < size; ++i)
    {
        const double magnitude = std::abs(values[i]);
        largest = std::max(largest, magnitude);
    }
    maximum = largest;
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
    const auto add_product = [&left, &right](std::size_t i, std::array<double, 1>& sum) {
        sum[0] += left[i] * right[i];
    };
    double value = thread_sums<1>(left.size(), add_product)[0];
    Parallel_Reduce::reduce_pool(value);
    return value * basis.omega / basis.nxyz;
}

// P r = eps^-1/2 G eps^-1/2 r; only FFT scratch survives an application.
void apply_precond(const std::vector<double>& rhs,
                   const std::vector<double>& invsqrt,
                   CoulombOperator& coulomb,
                   std::vector<double>& weighted,
                   std::vector<double>& value)
{
    const std::size_t size = rhs.size();
    weighted.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        weighted[i] = rhs[i] * invsqrt[i];
    }
    coulomb.apply_potential(weighted, value);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        value[i] *= invsqrt[i];
    }
}

// The preconditioned sqrt-CG system of one response.
struct SqrtCgOperator
{
    const std::vector<double>& coefficient;
    const std::vector<double>& invsqrt;
    const ModulePW::PW_Basis& basis;
    CoulombOperator& coulomb;
};

// Potential and charge residual of the sqrt-CG, with preconditioner scratch.
struct SqrtCgState
{
    std::vector<double> potential;
    std::vector<double> residual;
    std::vector<double> weighted;
    std::vector<double> preconditioned;
};

// Replace the cold start by one fixed-point update of the previous solution,
// phi = P (charge - c phi_old), when that lowers the RMS charge residual.
void try_warm_start(const SqrtCgOperator& system,
                    const std::vector<double>& charge,
                    const std::vector<double>& initial_potential,
                    SqrtCgState& state,
                    PolarizationResult& polarization)
{
    const std::size_t size = charge.size();
    const std::vector<double>& coefficient = system.coefficient;
    std::vector<double>& z = state.preconditioned;
    std::vector<double> guess_residual(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        guess_residual[i] = charge[i] - coefficient[i] * initial_potential[i];
    }
    apply_precond(guess_residual, system.invsqrt, system.coulomb, state.weighted, z);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        guess_residual[i] = coefficient[i] * (initial_potential[i] - z[i]);
    }
    double guess_rms = 0.0;
    double guess_max = 0.0;
    residual_norms(guess_residual, system.basis, guess_rms, guess_max);
    if (guess_rms < polarization.residual_rms)
    {
        state.potential.swap(z);
        state.residual.swap(guess_residual);
        polarization.warm_started = true;
        polarization.residual_rms = guess_rms;
        polarization.residual_max = guess_max;
    }
}

// Preconditioned CG until both residual tolerances hold; a breakdown or a
// missed tolerance stops the run with WARNING_QUIT.
void iterate_sqrt_cg(const SqrtCgOperator& system,
                     const PolarizationSolverParameters& solver,
                     SqrtCgState& state,
                     PolarizationResult& polarization)
{
    const ModulePW::PW_Basis& basis = system.basis;
    const std::vector<double>& coefficient = system.coefficient;
    std::vector<double>& potential = state.potential;
    std::vector<double>& residual = state.residual;
    std::vector<double>& z = state.preconditioned;
    const std::size_t size = residual.size();
    std::vector<double> direction(size, 0.0);
    std::vector<double> image(size, 0.0);
    double old_rz = 0.0;
    for (int iteration = 1; !converged(polarization, solver) && iteration <= solver.max_iterations; ++iteration)
    {
        apply_precond(residual, system.invsqrt, system.coulomb, state.weighted, z);
        const double rz = grid_dot(residual, z, basis);
        if (!std::isfinite(rz) || std::abs(rz) < 1e-30)
        {
            ModuleBase::WARNING_QUIT("ModuleSccs::solve_response",
                                     "SCCS sqrt-CG has a null or nonfinite preconditioned residual");
        }
        const double beta = std::abs(old_rz) > 1e-30 ? rz / old_rz : 0.0;
        old_rz = rz;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            direction[i] = z[i] + beta * direction[i];
            image[i] = coefficient[i] * z[i] + residual[i] + beta * image[i];
        }
        const double curvature = grid_dot(direction, image, basis);
        if (!std::isfinite(curvature) || curvature == 0.0)
        {
            ModuleBase::WARNING_QUIT("ModuleSccs::solve_response", "SCCS sqrt-CG has invalid curvature");
        }
        const double alpha = rz / curvature;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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
        ModuleBase::WARNING_QUIT("ModuleSccs::solve_response", message.str());
    }
}

// Keep the unshifted solution for warm starts, fix the gauge of a periodic
// cell and store the potential, its gradient and the boundary potentials.
void finish_response(const SqrtCgOperator& system,
                     const std::vector<double>& charge,
                     const CavityParameters& cavity,
                     const CavityDerivatives& derivatives,
                     double tpiba,
                     std::vector<double>& potential,
                     SccsResponse& response)
{
    const ModulePW::PW_Basis& basis = system.basis;
    const std::vector<double>& coefficient = system.coefficient;
    const std::vector<double>& invsqrt = system.invsqrt;
    const bool open_boundary = system.coulomb.has_correction();
    const std::size_t size = potential.size();
    PolarizationResult& polarization = response.polarization;
    response.restart_potential = potential;
    // Preserve the physical gauge fixed by a boundary correction.
    if (!open_boundary)
    {
        const auto add_value = [&potential](std::size_t i, std::array<double, 1>& sum) {
            sum[0] += potential[i];
        };
        double mean = thread_sums<1>(size, add_value)[0];
        Parallel_Reduce::reduce_pool(mean);
        mean /= basis.nxyz;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            potential[i] -= mean;
        }
    }
    // The same condition under which make_cav_derivs keeps the filtered gradient.
    const bool lowpass = open_boundary && uses_lowpass(cavity);
    if (!lowpass)
    {
        std::vector<std::complex<double>> potential_g(basis.npw);
        basis.real2recip(potential.data(), potential_g.data());
        polarization.gradient.resize(size);
        XC_Functional::grad_rho(potential_g.data(), polarization.gradient.data(), &basis, tpiba);
    }
    polarization.potential.swap(potential);
    if (open_boundary)
    {
        finish_open_bc(charge, coefficient, invsqrt, cavity, basis, response);
    }
    if (lowpass)
    {
        lowpass_bnd_pot(charge, cavity, derivatives, basis, tpiba, response);
    }
    else
    {
        // Continuum derivative -(d eps/ds) |grad v|^2/(8 pi), as in Environ: exact
        // for the discrete energy only when |grad v|^2 is resolved by the FFT grid.
        const double log_bulk = std::log(cavity.epsilon_bulk);
        response.boundary_potential.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            const ModuleBase::Vector3<double>& gradient = polarization.gradient[i];
            const double gradient_square = gradient * gradient;
            response.boundary_potential[i]
                = log_bulk * response.epsilon[i] * gradient_square / (8.0 * ModuleBase::PI);
        }
    }
}
} // namespace

void boundary_to_rho(const SccsResponse& response,
                     const std::vector<double>& probe_kernel,
                     const ModulePW::PW_Basis& basis,
                     const std::vector<double>& boundary_potential,
                     std::vector<double>& density_potential)
{
    // Environ calculator: sa_de_dboundary, then de_dboundary * dscaled.
    if (response.solvent_aware.enabled())
    {
        density_potential = sa_adjoint(response.solvent_aware, probe_kernel, basis, boundary_potential);
    }
    else
    {
        density_potential = boundary_potential;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < density_potential.size(); ++i)
    {
        density_potential[i] *= response.dsolute_drho[i];
    }
}

void solve_response(const std::vector<double>& density,
                    const std::vector<double>& charge,
                    const CavityParameters& cavity,
                    const PolarizationSolverParameters& solver,
                    const std::vector<double>& initial_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    SccsResponse& result)
{
    PeriodicCoulombOperator coulomb(basis, tpiba);
    const std::vector<double> no_probe;
    solve_response(density, charge, cavity, no_probe, solver, initial_potential, basis, tpiba, coulomb, result);
}

void solve_response(const std::vector<double>& density,
                    const std::vector<double>& charge,
                    const CavityParameters& cavity,
                    const std::vector<double>& probe_kernel,
                    const PolarizationSolverParameters& solver,
                    const std::vector<double>& initial_potential,
                    const ModulePW::PW_Basis& basis,
                    double tpiba,
                    CoulombOperator& coulomb,
                    SccsResponse& result)
{
    // Determine start mode collectively, including ranks with no real-space points.
    double initial_count = initial_potential.size();
    Parallel_Reduce::reduce_pool(initial_count);
    const bool warm_start = initial_count != 0.0;

    const bool open_boundary = coulomb.has_correction();
    SccsResponse candidate;
    CavityDerivatives derivatives;
    make_cav_derivs(density, cavity, probe_kernel, basis, tpiba, open_boundary, candidate, derivatives);
    const std::vector<double>& coefficient = derivatives.coefficient;
    const std::size_t size = density.size();
    std::vector<double> invsqrt(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        invsqrt[i] = 1.0 / std::sqrt(candidate.epsilon[i]);
    }
    const SqrtCgOperator system = {coefficient, invsqrt, basis, coulomb};
    SqrtCgState state;
    state.residual = charge;
    state.potential.assign(size, 0.0);
    PolarizationResult& polarization = candidate.polarization;
    residual_norms(state.residual, basis, polarization.residual_rms, polarization.residual_max);
    if (!converged(polarization, solver) && warm_start)
    {
        try_warm_start(system, charge, initial_potential, state, polarization);
    }
    iterate_sqrt_cg(system, solver, state, polarization);
    std::vector<double>& potential = state.potential;
    finish_response(system, charge, cavity, derivatives, tpiba, potential, candidate);
    boundary_to_rho(candidate, probe_kernel, basis, candidate.boundary_potential,
                    candidate.cavity_potential);
    result = std::move(candidate);
}
} // namespace ModuleSccs
