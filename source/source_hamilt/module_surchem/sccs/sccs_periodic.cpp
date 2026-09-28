#include "sccs_periodic.h"

#include "sccs_pw_coulomb.h"
#include "sccs_pw_reduction.h"

#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace ModuleSccs
{

PeriodicSccsResult solve_periodic_sccs(
    const std::vector<double>& cavity_density,
    const std::vector<double>& solute_charge,
    const CavityParameters& cavity_parameters,
    const PolarizationSolverParameters& solver_parameters,
    const std::vector<double>& initial_polarization_charge,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const int pool_process_count)
{
    const PeriodicCoulombOperator coulomb(basis, tpiba);
    const PoolPolarizationReduction reduction(pool_process_count);
    return solve_sccs_response(cavity_density,
                               solute_charge,
                               cavity_parameters,
                               solver_parameters,
                               initial_polarization_charge,
                               basis,
                               tpiba,
                               coulomb,
                               reduction);
}

PeriodicSccsResult solve_sccs_response(
    const std::vector<double>& density,
    const std::vector<double>& charge,
    const ModuleSccs::CavityParameters& cavity,
    const ModuleSccs::PolarizationSolverParameters& solver,
    const std::vector<double>& initial,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const ModuleSccs::CoulombOperator& coulomb,
    const ModuleSccs::PolarizationReduction& reduction)
{
    ModuleSccs::PeriodicSccsResult result;
    const auto density_gradient = ModuleSccs::periodic_gradient(density, basis, tpiba);
    result.density_gradient = density_gradient;
    const std::size_t size = density.size();
    result.solute.resize(size);
    result.dsolute_drho.resize(size);
    result.epsilon.resize(size);
    result.depsilon_drho.resize(size);
    result.grad_log_epsilon.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const auto point = ModuleSccs::evaluate_cavity(density[i], cavity);
        result.solute[i] = point.solute;
        result.dsolute_drho[i] = point.dsolute_drho;
        result.epsilon[i] = point.epsilon;
        result.depsilon_drho[i] = point.depsilon_drho;
        const double coefficient = point.depsilon_drho / point.epsilon;
        for (int d = 0; d < 3; ++d)
            result.grad_log_epsilon[i][d] = coefficient * density_gradient[i][d];
    }
    // Environ dielectric::factsqrt for electronic chain derivatives, in Ha units.
    std::vector<std::complex<double>> density_g(basis.npw);
    basis.real2recip(density.data(), density_g.data());
    for (int ig = 0; ig < basis.npw; ++ig)
        density_g[ig] *= -tpiba * tpiba * basis.gg[ig];
    std::vector<double> laplacian(size);
    basis.recip2real(density_g.data(), laplacian.data());
    std::vector<double> coefficient(size);
    std::vector<double> invsqrt(size);
    const double density_ratio = cavity.density_max / cavity.density_min;
    const double width = std::log(density_ratio);
    const double log_bulk = std::log(cavity.epsilon_bulk);
    for (std::size_t i = 0; i < size; ++i)
    {
        double second_log = 0.0;
        if (density[i] > cavity.density_min && density[i] < cavity.density_max)
        {
            const double local_density_ratio = cavity.density_max / density[i];
            const double x = std::log(local_density_ratio)/width;
            const double angle = ModuleBase::TWO_PI*x;
            second_log = log_bulk*(1.0-std::cos(angle)
                         + ModuleBase::TWO_PI*std::sin(angle)/width)
                         /(width*density[i]*density[i]);
        }
        const double first_log = result.depsilon_drho[i]/result.epsilon[i];
        double gradient_square = 0.0;
        for (int d = 0; d < 3; ++d)
            gradient_square += density_gradient[i][d]*density_gradient[i][d];
        const double lap_log = first_log*laplacian[i]+second_log*gradient_square;
        coefficient[i] = result.epsilon[i]*(0.5*lap_log
                         + 0.25*first_log*first_log*gradient_square)/ModuleBase::FOUR_PI;
        invsqrt[i] = 1.0/std::sqrt(result.epsilon[i]);
    }
    // P r = epsilon^-1/2 C_PCC epsilon^-1/2 r.
    // Per-solve scratch is shared by CG applications; the cavity is fixed
    // throughout this solve. Only the scalar Coulomb potential is needed.
    std::vector<double> weighted(size);
    std::vector<double> potential_work;
    const auto precondition = [&](const std::vector<double>& rhs, std::vector<double>& value)
    {
        for (std::size_t i = 0; i < size; ++i) weighted[i] = rhs[i]*invsqrt[i];
        coulomb.apply_potential(weighted, potential_work);
        value.resize(size);
        for (std::size_t i = 0; i < size; ++i) value[i] = potential_work[i]*invsqrt[i];
    };
    const auto dot = [&](const std::vector<double>& left, const std::vector<double>& right)
    {
        double value = 0.0;
        for (std::size_t i = 0; i < size; ++i) value += left[i]*right[i];
        reduction.reduce_sum(value);
        return value * basis.omega / basis.nxyz;
    };
    std::vector<double> residual = charge;
    std::vector<double> potential(size, 0.0);
    std::vector<double> direction(size, 0.0);
    std::vector<double> image(size, 0.0);
    std::vector<double> z;
    double old_rz = 0.0;
    double initial_square = 0.0;
    for (const double value : residual)
    {
        initial_square += value * value;
    }
    reduction.reduce_sum(initial_square);
    bool converged = initial_square <= solver.tolerance_rms;
    for (int iteration = 1; !converged && iteration <= solver.max_iterations; ++iteration)
    {
        precondition(residual, z);
        const double rz = dot(residual, z);
        if (!std::isfinite(rz) || std::abs(rz) < 1e-30)
            throw std::runtime_error("CG sqrt null/nonfinite preconditioned residual");
        const double beta = std::abs(old_rz) > 1e-30 ? rz/old_rz : 0.0;
        old_rz = rz;
        for (std::size_t i = 0; i < size; ++i)
        {
            direction[i] = z[i]+beta*direction[i];
            image[i] = coefficient[i]*z[i]+residual[i]+beta*image[i];
        }
        const double curvature = dot(direction, image);
        if (!std::isfinite(curvature) || curvature == 0.0)
            throw std::runtime_error("CG sqrt invalid curvature");
        const double alpha = rz/curvature;
        double square = 0.0;
        double maximum = 0.0;
        double count = size;
        for (std::size_t i = 0; i < size; ++i)
        {
            potential[i] += alpha*direction[i];
            residual[i] -= alpha*image[i];
            square += residual[i]*residual[i];
            maximum = std::max(maximum, std::abs(residual[i]));
        }
        reduction.reduce_residual(square, maximum, count);
        result.polarization.iterations = iteration;
        const double mean_square = square / count;
        result.polarization.residual_rms = std::sqrt(mean_square);
        result.polarization.residual_max = maximum;
        if (square <= solver.tolerance_rms)
        {
            converged = true;
            break;
        }
        if (!std::isfinite(square)) throw std::runtime_error("CG sqrt diverged");
    }
    if (!converged) throw std::runtime_error("Experimental CG sqrt reached iteration limit");
    // Independently check the preconditioned equation v=P(q-Kv).
    std::vector<double> right(size);
    for (std::size_t i = 0; i < size; ++i) right[i] = charge[i]-coefficient[i]*potential[i];
    precondition(right, z);
    double square = 0.0;
    double maximum = 0.0;
    double count = size;
    for (std::size_t i = 0; i < size; ++i)
    {
        const double defect = potential[i]-z[i];
        square += defect*defect;
        maximum = std::max(maximum, std::abs(defect));
    }
    reduction.reduce_residual(square, maximum, count);
    const double defect_mean_square = square / count;
    if (basis.startz_current == 0)
        std::cout << "CG_FIXED_POINT_DEFECT " << std::sqrt(defect_mean_square)
                  << " " << maximum << std::endl;
    // Recover induced charge for existing diagnostic/state consumers.
    basis.real2recip(potential.data(), density_g.data());
    for (int ig = 0; ig < basis.npw; ++ig)
        density_g[ig] *= tpiba*tpiba*basis.gg[ig]/ModuleBase::FOUR_PI;
    result.polarization.polarization_charge.resize(size);
    basis.recip2real(density_g.data(), result.polarization.polarization_charge.data());
    for (std::size_t i = 0; i < size; ++i)
        result.polarization.polarization_charge[i] -= charge[i];
    result.polarization.status = ModuleSccs::PolarizationStatus::Converged;
    double mean = 0.0;
    for (double value : potential) mean += value;
    reduction.reduce_sum(mean);
    mean /= basis.nxyz;
    for (double& value : potential) value -= mean;
    result.polarization.field.potential = potential;
    // Environ dielectric::de_dboundary differentiates the solved potential on
    // its derivative grid, including when the potential contains a PCC term.
    result.polarization.field.gradient = ModuleSccs::periodic_gradient(potential, basis, tpiba);
    return result;
}

} // namespace ModuleSccs
