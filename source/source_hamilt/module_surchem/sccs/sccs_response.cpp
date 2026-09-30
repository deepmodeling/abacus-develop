#include "sccs_response.h"

#include "sccs_pw_coulomb.h"
#include "../common/charge_reduction.h"

#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace ModuleSccs
{

std::vector<double> continuum_polarization_charge(
    const std::vector<double>& solute_charge,
    const SccsResponse& response)
{
    const std::size_t size = solute_charge.size();
    if (response.epsilon.size() != size || response.grad_log_epsilon.size() != size
        || response.polarization.field.gradient.size() != size)
    {
        throw std::invalid_argument("SCCS polarization source arrays must have the same size");
    }
    std::vector<double> polarization(size);
    for (std::size_t index = 0; index < size; ++index)
    {
        const double epsilon = response.epsilon[index];
        if (!std::isfinite(epsilon) || epsilon < 1.0)
        {
            throw std::domain_error("SCCS polarization source requires finite epsilon >= 1");
        }
        double projection = 0.0;
        for (int direction = 0; direction < 3; ++direction)
        {
            projection += response.grad_log_epsilon[index][direction]
                          * response.polarization.field.gradient[index][direction];
        }
        polarization[index] = projection / ModuleBase::FOUR_PI
                              + solute_charge[index] * (1.0 / epsilon - 1.0);
        if (!std::isfinite(polarization[index]))
        {
            throw std::domain_error("SCCS polarization source must be finite");
        }
    }
    return polarization;
}

namespace
{

SccsResponse prepare_cavity(
    const std::vector<double>& density,
    const CavityParameters& cavity,
    const ModulePW::PW_Basis& basis,
    const double tpiba)
{
    SccsResponse result;
    const std::vector<ModuleBase::Vector3<double>> density_gradient
        = ModuleSccs::periodic_gradient(density, basis, tpiba);
    result.density_gradient = density_gradient;
    const std::size_t size = density.size();
    result.solute.resize(size);
    result.dsolute_drho.resize(size);
    result.epsilon.resize(size);
    result.depsilon_drho.resize(size);
    result.grad_log_epsilon.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const CavityPoint point = ModuleSccs::evaluate_cavity(density[i], cavity);
        result.solute[i] = point.solute;
        result.dsolute_drho[i] = point.dsolute_drho;
        result.epsilon[i] = point.epsilon;
        result.depsilon_drho[i] = point.depsilon_drho;
        const double coefficient = point.depsilon_drho / point.epsilon;
        for (int d = 0; d < 3; ++d)
        {
            result.grad_log_epsilon[i][d] = coefficient * density_gradient[i][d];
        }
    }
    return result;
}

// Pool RMS and maximum absolute value of a distributed grid array.
void reduced_rms_max(const std::vector<double>& values,
                     const ModuleSurchem::ChargeReduction& reduction,
                     double& rms,
                     double& maximum)
{
    double square = 0.0;
    double local_maximum = 0.0;
    double count = static_cast<double>(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        square += values[i] * values[i];
        local_maximum = std::max(local_maximum, std::abs(values[i]));
    }
    reduction.reduce_sum(square);
    reduction.reduce_max(local_maximum);
    reduction.reduce_sum(count);
    if (!std::isfinite(square) || !std::isfinite(local_maximum) || count <= 0.0)
    {
        throw std::runtime_error("SCCS sqrt-CG residual is not finite");
    }
    rms = std::sqrt(square / count);
    maximum = local_maximum;
}

// Update the residual norms of polarization; true when both pass
// sccs_tol_rms and sccs_tol_max.
bool residual_converged(const std::vector<double>& residual,
                        const PolarizationSolverParameters& solver,
                        const ModuleSurchem::ChargeReduction& reduction,
                        PolarizationResult& polarization)
{
    reduced_rms_max(residual, reduction, polarization.residual_rms, polarization.residual_max);
    return polarization.residual_rms <= solver.tolerance_rms
           && polarization.residual_max <= solver.tolerance_max;
}

// Grid inner product of two distributed arrays.
double grid_dot(const std::vector<double>& left,
                const std::vector<double>& right,
                const ModulePW::PW_Basis& basis,
                const ModuleSurchem::ChargeReduction& reduction)
{
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        value += left[i] * right[i];
    }
    reduction.reduce_sum(value);
    return value * basis.omega / basis.nxyz;
}

// P r = eps^-1/2 G eps^-1/2 r with the periodic or PCC-corrected Coulomb
// operator G. The cavity is fixed during a solve, so the scratch arrays are
// reused by every application; only the scalar potential is needed.
class SqrtPreconditioner
{
  public:
    SqrtPreconditioner(const std::vector<double>& invsqrt, const CoulombOperator& coulomb)
        : invsqrt_(invsqrt), coulomb_(coulomb)
    {
        this->weighted_.resize(invsqrt.size());
    }

    void apply(const std::vector<double>& rhs, std::vector<double>& value)
    {
        const std::size_t size = this->invsqrt_.size();
        for (std::size_t i = 0; i < size; ++i)
        {
            this->weighted_[i] = rhs[i] * this->invsqrt_[i];
        }
        this->coulomb_.apply_potential(this->weighted_, this->potential_);
        value.resize(size);
        for (std::size_t i = 0; i < size; ++i)
        {
            value[i] = this->potential_[i] * this->invsqrt_[i];
        }
    }

  private:
    const std::vector<double>& invsqrt_;
    const CoulombOperator& coulomb_;
    std::vector<double> weighted_;
    std::vector<double> potential_;
};

// Environ dielectric::factsqrt for electronic chain derivatives, in Ha units.
void chain_factsqrt(const std::vector<double>& density,
                    const CavityParameters& cavity,
                    const SccsResponse& result,
                    const ModulePW::PW_Basis& basis,
                    const double tpiba,
                    std::vector<double>& coefficient)
{
    const std::size_t size = density.size();
    const std::vector<ModuleBase::Vector3<double>>& density_gradient = result.density_gradient;
    std::vector<std::complex<double>> density_g(basis.npw);
    basis.real2recip(density.data(), density_g.data());
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        density_g[ig] *= -tpiba * tpiba * basis.gg[ig];
    }
    std::vector<double> laplacian(size);
    basis.recip2real(density_g.data(), laplacian.data());
    const double density_ratio = cavity.density_max / cavity.density_min;
    const double width = std::log(density_ratio);
    const double log_bulk = std::log(cavity.epsilon_bulk);
    for (std::size_t i = 0; i < size; ++i)
    {
        double second_log = 0.0;
        if (density[i] > cavity.density_min && density[i] < cavity.density_max)
        {
            const double local_density_ratio = cavity.density_max / density[i];
            const double x = std::log(local_density_ratio) / width;
            const double angle = ModuleBase::TWO_PI * x;
            second_log = log_bulk * (1.0 - std::cos(angle) + ModuleBase::TWO_PI * std::sin(angle) / width)
                         / (width * density[i] * density[i]);
        }
        const double first_log = result.depsilon_drho[i] / result.epsilon[i];
        double gradient_square = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            gradient_square += density_gradient[i][d] * density_gradient[i][d];
        }
        const double lap_log = first_log * laplacian[i] + second_log * gradient_square;
        coefficient[i] = result.epsilon[i] * (0.5 * lap_log + 0.25 * first_log * first_log * gradient_square)
                         / ModuleBase::FOUR_PI;
    }
}

// Environ 3.1.1 core_fft_lowpass: with lowpass_p1 and lowpass_p2 positive,
// every switching-function derivative is multiplied by
// 0.5 erfc(p1 G^2/Gcut^2 - p2), Gcut^2 being the density cutoff; else by one.
double switching_filter(const double gg,
                        const CavityParameters& cavity,
                        const ModulePW::PW_Basis& basis)
{
    if (!uses_switching_lowpass(cavity))
    {
        return 1.0;
    }
    const double argument = cavity.lowpass_p1 * gg / basis.ggecut - cavity.lowpass_p2;
    return 0.5 * std::erfc(argument);
}

// Spectral gradient of the switching function, filtered when requested.
std::vector<ModuleBase::Vector3<double>> switching_gradient(const std::vector<double>& values,
                                                            const CavityParameters& cavity,
                                                            const ModulePW::PW_Basis& basis,
                                                            const double tpiba)
{
    if (!uses_switching_lowpass(cavity))
    {
        return ModuleSccs::periodic_gradient(values, basis, tpiba);
    }
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    std::vector<std::complex<double>> gradient_g(basis.npw);
    std::vector<double> gradient_r(values.size());
    std::vector<ModuleBase::Vector3<double>> gradient(values.size());
    for (int d = 0; d < 3; ++d)
    {
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double filter = switching_filter(basis.gg[ig], cavity, basis);
            gradient_g[ig] = ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][d] * values_g[ig] * filter;
        }
        basis.recip2real(gradient_g.data(), gradient_r.data());
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            gradient[i][d] = gradient_r[i];
        }
    }
    return gradient;
}

// Spectral Laplacian of the switching function. It is symmetric under the
// grid inner product, so it is its own transpose in the cavity derivative.
void switching_laplacian(const std::vector<double>& values,
                         const CavityParameters& cavity,
                         const ModulePW::PW_Basis& basis,
                         const double tpiba,
                         std::vector<double>& laplacian)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const double filter = switching_filter(basis.gg[ig], cavity, basis);
        values_g[ig] *= -tpiba * tpiba * basis.gg[ig] * filter;
    }
    laplacian.resize(values.size());
    basis.recip2real(values_g.data(), laplacian.data());
}

// Spectral divergence matching switching_gradient; minus this divergence is
// the transpose of switching_gradient.
void switching_divergence(const std::vector<ModuleBase::Vector3<double>>& field,
                          const CavityParameters& cavity,
                          const ModulePW::PW_Basis& basis,
                          const double tpiba,
                          std::vector<double>& divergence)
{
    const std::size_t size = field.size();
    std::vector<double> component(size);
    std::vector<std::complex<double>> component_g(basis.npw);
    const std::complex<double> zero(0.0, 0.0);
    std::vector<std::complex<double>> divergence_g(basis.npw, zero);
    for (int d = 0; d < 3; ++d)
    {
        for (std::size_t i = 0; i < size; ++i)
        {
            component[i] = field[i][d];
        }
        basis.real2recip(component.data(), component_g.data());
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            const double filter = switching_filter(basis.gg[ig], cavity, basis);
            divergence_g[ig] += ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][d] * component_g[ig] * filter;
        }
    }
    divergence.resize(size);
    basis.recip2real(divergence_g.data(), divergence.data());
}

// With PCC the potential at the cavity edge carries the open-boundary monopole
// and dipole. The chain-rule f drops steeply to zero at density_min, and the
// sampled f v source then fails to converge with the grid. Differentiate the
// switching function s on the FFT grid instead (Environ deriv_method 'fft'):
// ln eps = ln(eps_bulk) (1 - s). grad ln(eps) is replaced consistently, and
// grad s is returned for the cavity derivative.
void switching_fft_factsqrt(const CavityParameters& cavity,
                            const ModulePW::PW_Basis& basis,
                            const double tpiba,
                            SccsResponse& result,
                            std::vector<double>& coefficient,
                            std::vector<ModuleBase::Vector3<double>>& solute_gradient)
{
    const std::size_t size = result.solute.size();
    const double log_bulk = std::log(cavity.epsilon_bulk);
    solute_gradient = switching_gradient(result.solute, cavity, basis, tpiba);
    std::vector<double> solute_laplacian;
    switching_laplacian(result.solute, cavity, basis, tpiba, solute_laplacian);
    for (std::size_t i = 0; i < size; ++i)
    {
        double gradient_square = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            result.grad_log_epsilon[i][d] = -log_bulk * solute_gradient[i][d];
            gradient_square += solute_gradient[i][d] * solute_gradient[i][d];
        }
        const double lap_log = -log_bulk * solute_laplacian[i];
        coefficient[i] = result.epsilon[i]
                         * (0.5 * lap_log + 0.25 * log_bulk * log_bulk * gradient_square)
                         / ModuleBase::FOUR_PI;
    }
}

// Continuum cavity potential -eps'|grad v|^2/(8 pi) of Environ
// dielectric::de_dboundary, with grad v from the solved potential.
void continuum_cavity_potential(SccsResponse& result)
{
    const std::size_t size = result.depsilon_drho.size();
    result.cavity_potential.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const ModuleBase::Vector3<double>& gradient = result.polarization.field.gradient[i];
        const double gradient_square
            = gradient.x * gradient.x + gradient.y * gradient.y + gradient.z * gradient.z;
        result.cavity_potential[i]
            = -(result.depsilon_drho[i] * gradient_square / (8.0 * ModuleBase::PI));
    }
}

// Exact derivative of the discrete reaction energy through the cavity for the
// lowpass switching-function factsqrt. The sqrt-CG solves A v = q with
// A = sqrt(eps) G^-1 sqrt(eps) + F, which is symmetric, so E = q^T A^-1 q / 2
// changes by dE = -v^T dA v / 2 and needs no adjoint solve. With
// L = ln(eps_bulk), s the switching function and b = eps v^2/(8 pi), the
// sqrt(eps) term and the pointwise part of dF add up to q v / 2, and the
// transposes of the filtered lapl and grad in F give
// dE/ds = L/2 (q v + lapl b + L div(b grad s)).
// Its continuum limit is -eps'|grad v|^2/(8 pi). Without the filter the
// discrete derivative is not grid-converged at the cavity edge and its
// grid-scale oscillations drive the electronic SCF to diverge.
void switching_cavity_potential(const std::vector<double>& charge,
                                const std::vector<double>& potential,
                                const std::vector<ModuleBase::Vector3<double>>& solute_gradient,
                                const CavityParameters& cavity,
                                const ModulePW::PW_Basis& basis,
                                const double tpiba,
                                SccsResponse& result)
{
    const std::size_t size = charge.size();
    const double log_bulk = std::log(cavity.epsilon_bulk);
    std::vector<double> weight(size);
    std::vector<ModuleBase::Vector3<double>> weighted_gradient(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        weight[i] = result.epsilon[i] * potential[i] * potential[i] / (8.0 * ModuleBase::PI);
        for (int d = 0; d < 3; ++d)
        {
            weighted_gradient[i][d] = weight[i] * solute_gradient[i][d];
        }
    }
    std::vector<double> weight_laplacian;
    switching_laplacian(weight, cavity, basis, tpiba, weight_laplacian);
    std::vector<double> weighted_divergence;
    switching_divergence(weighted_gradient, cavity, basis, tpiba, weighted_divergence);
    result.cavity_potential.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const double derivative = 0.5 * log_bulk
                                  * (charge[i] * potential[i] + weight_laplacian[i]
                                     + log_bulk * weighted_divergence[i]);
        result.cavity_potential[i] = derivative * result.dsolute_drho[i];
    }
}

// Check the preconditioned equation v = P(q - K v) independently of the CG
// recurrences; this costs one extra Poisson solve.
void check_fixed_point(const std::vector<double>& charge,
                       const std::vector<double>& coefficient,
                       const std::vector<double>& potential,
                       const ModuleSurchem::ChargeReduction& reduction,
                       SqrtPreconditioner& preconditioner,
                       PolarizationResult& polarization)
{
    const std::size_t size = potential.size();
    std::vector<double> right(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        right[i] = charge[i] - coefficient[i] * potential[i];
    }
    std::vector<double> image;
    preconditioner.apply(right, image);
    std::vector<double> defect(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        defect[i] = potential[i] - image[i];
    }
    reduced_rms_max(defect, reduction,
                    polarization.fixed_point_defect_rms,
                    polarization.fixed_point_defect_max);
    polarization.fixed_point_checked = true;
}

// Shift a periodic potential to zero cell mean (ENVIRON generalized_sqrt).
void remove_mean(const ModulePW::PW_Basis& basis,
                 const ModuleSurchem::ChargeReduction& reduction,
                 std::vector<double>& potential)
{
    double mean = 0.0;
    for (const double value : potential)
    {
        mean += value;
    }
    reduction.reduce_sum(mean);
    mean /= basis.nxyz;
    for (double& value : potential)
    {
        value -= mean;
    }
}

// The CG builds sqrt(eps) v = w = C_PCC(s) with s = (q - f v)/sqrt(eps).
// Report the ENVIRON dielectric_of_potential polarization density and the
// far-field polarization charge int(s)/sqrt(eps_bulk) - int(q).
void finish_open_boundary_response(const std::vector<double>& charge,
                                   const std::vector<double>& coefficient,
                                   const std::vector<double>& invsqrt,
                                   const CavityParameters& cavity,
                                   const ModulePW::PW_Basis& basis,
                                   const ModuleSurchem::ChargeReduction& reduction,
                                   SccsResponse& result)
{
    const std::size_t size = charge.size();
    const std::vector<double>& potential = result.polarization.field.potential;
    double source_sum = 0.0;
    double solute_sum = 0.0;
    for (std::size_t i = 0; i < size; ++i)
    {
        source_sum += (charge[i] - coefficient[i] * potential[i]) * invsqrt[i];
        solute_sum += charge[i];
    }
    // The corrected potential has no periodic Laplacian inverse; use the
    // ENVIRON dielectric_of_potential polarization charge instead. Its integral
    // carries a finite-grid error; the far field fixes the net screening charge.
    result.polarization.polarization_charge = continuum_polarization_charge(charge, result);
    reduction.reduce_sum(source_sum);
    reduction.reduce_sum(solute_sum);
    const double volume_element = basis.omega / basis.nxyz;
    const double bulk_invsqrt = 1.0 / std::sqrt(cavity.epsilon_bulk);
    result.far_field_polarization_charge
        = (source_sum * bulk_invsqrt - solute_sum) * volume_element;
}

} // namespace

SccsResponse solve_sccs_response(
    const std::vector<double>& density,
    const std::vector<double>& charge,
    const ModuleSccs::CavityParameters& cavity,
    const ModuleSccs::PolarizationSolverParameters& solver,
    const std::vector<double>& initial_potential,
    const ModulePW::PW_Basis& basis,
    const double tpiba,
    const ModuleSccs::CoulombOperator& coulomb,
    const ModuleSurchem::ChargeReduction& reduction)
{
    ModuleBase::timer::start("ModuleSccs", "solve_sccs_response");
    ModuleSccs::validate_cavity_parameters(cavity);
    SccsResponse result = prepare_cavity(density, cavity, basis, tpiba);
    const std::size_t size = density.size();
    std::vector<double> coefficient(size);
    std::vector<double> invsqrt(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        invsqrt[i] = 1.0 / std::sqrt(result.epsilon[i]);
    }
    const bool open_boundary = coulomb.has_boundary_correction();
    std::vector<ModuleBase::Vector3<double>> solute_gradient;
    if (!open_boundary)
    {
        chain_factsqrt(density, cavity, result, basis, tpiba, coefficient);
    }
    else
    {
        switching_fft_factsqrt(cavity, basis, tpiba, result, coefficient, solute_gradient);
    }
    SqrtPreconditioner preconditioner(invsqrt, coulomb);
    std::vector<double> residual = charge;
    std::vector<double> potential(size, 0.0);
    std::vector<double> direction(size, 0.0);
    std::vector<double> image(size, 0.0);
    std::vector<double> z;
    double old_rz = 0.0;
    PolarizationResult& polarization = result.polarization;
    bool converged = residual_converged(residual, solver, reduction, polarization);
    // ENVIRON generalized_sqrt warm start: one preconditioned fixed-point step
    // v = P(q - K v_old) from the previous potential, whose charge residual is
    // K (v_old - v). Keep it only when it improves on the cold-start residual.
    if (!converged && initial_potential.size() == size)
    {
        std::vector<double> guess_residual(size);
        for (std::size_t i = 0; i < size; ++i)
        {
            guess_residual[i] = charge[i] - coefficient[i] * initial_potential[i];
        }
        preconditioner.apply(guess_residual, z);
        for (std::size_t i = 0; i < size; ++i)
        {
            guess_residual[i] = coefficient[i] * (initial_potential[i] - z[i]);
        }
        double guess_rms = 0.0;
        double guess_max = 0.0;
        reduced_rms_max(guess_residual, reduction, guess_rms, guess_max);
        if (guess_rms < polarization.residual_rms)
        {
            potential.swap(z);
            residual.swap(guess_residual);
            polarization.warm_started = true;
            converged = residual_converged(residual, solver, reduction, polarization);
        }
    }
    for (int iteration = 1; !converged && iteration <= solver.max_iterations; ++iteration)
    {
        preconditioner.apply(residual, z);
        const double rz = grid_dot(residual, z, basis, reduction);
        if (!std::isfinite(rz) || std::abs(rz) < 1e-30)
        {
            throw std::runtime_error("CG sqrt null/nonfinite preconditioned residual");
        }
        const double beta = std::abs(old_rz) > 1e-30 ? rz / old_rz : 0.0;
        old_rz = rz;
        for (std::size_t i = 0; i < size; ++i)
        {
            direction[i] = z[i] + beta * direction[i];
            image[i] = coefficient[i] * z[i] + residual[i] + beta * image[i];
        }
        const double curvature = grid_dot(direction, image, basis, reduction);
        if (!std::isfinite(curvature) || curvature == 0.0)
        {
            throw std::runtime_error("CG sqrt invalid curvature");
        }
        const double alpha = rz / curvature;
        for (std::size_t i = 0; i < size; ++i)
        {
            potential[i] += alpha * direction[i];
            residual[i] -= alpha * image[i];
        }
        polarization.iterations = iteration;
        converged = residual_converged(residual, solver, reduction, polarization);
    }
    if (!converged)
    {
        throw std::runtime_error(
            "SCCS sqrt-CG did not reach sccs_tol_rms and sccs_tol_max within sccs_maxiter");
    }
    if (solver.check_fixed_point)
    {
        check_fixed_point(charge, coefficient, potential, reduction, preconditioner, polarization);
    }
    result.restart_potential = potential;
    // A PCC operator in the preconditioner fixes the physical gauge; only the
    // periodic potential is shifted to zero mean.
    if (!open_boundary)
    {
        remove_mean(basis, reduction, potential);
    }
    result.polarization.field.potential = potential;
    // Environ dielectric::de_dboundary differentiates the solved potential on
    // its derivative grid, for the continuum cavity potential and diagnostics.
    result.polarization.field.gradient = ModuleSccs::periodic_gradient(potential, basis, tpiba);
    if (open_boundary)
    {
        finish_open_boundary_response(charge, coefficient, invsqrt, cavity, basis, reduction,
                                      result);
    }
    if (uses_switching_lowpass(cavity))
    {
        switching_cavity_potential(charge, potential, solute_gradient, cavity, basis, tpiba,
                                   result);
    }
    else
    {
        continuum_cavity_potential(result);
    }
    ModuleBase::timer::end("ModuleSccs", "solve_sccs_response");
    return result;
}

} // namespace ModuleSccs
