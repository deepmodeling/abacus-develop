#include "sccs_cav_deriv.h"
#include "sccs_cavity.h"
#include "sccs_lowpass.h"
#include "sccs_response.h"
#include "sccs_sa.h"
#include "sccs_thread_sum.h"

#include "source_base/constants.h"
#include "source_base/parallel_reduce.h"
#include "source_base/tool_quit.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
namespace
{
// Gradient and Laplacian of the field with coefficients values_g.
void spec_derivs(const std::vector<std::complex<double>>& values_g,
                 const ModulePW::PW_Basis& basis,
                 double tpiba,
                 std::vector<ModuleBase::Vector3<double>>& gradient,
                 std::vector<double>& laplacian)
{
    gradient.resize(basis.nrxx);
    laplacian.resize(basis.nrxx);
    XC_Functional::grad_rho(values_g.data(), gradient.data(), &basis, tpiba);
    XC_Functional::laplacian_rho(values_g.data(), laplacian.data(), &basis, tpiba);
}

// Gradient and Laplacian on the FFT grid, multiplied by filter unless it is empty.
void fft_derivatives(const std::vector<double>& values,
                     const std::vector<double>& filter,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     std::vector<ModuleBase::Vector3<double>>& gradient,
                     std::vector<double>& laplacian)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    if (!filter.empty())
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            values_g[ig] *= filter[ig];
        }
    }
    spec_derivs(values_g, basis, tpiba, gradient, laplacian);
}

// FFT derivatives of the local s and of the probe solute fraction c = p * s from
// one transform of s: the probe is diagonal in G, so p * grad s = grad c.
void fft_fill_derivs(const std::vector<double>& probe_kernel,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     const SccsResponse& response,
                     std::vector<ModuleBase::Vector3<double>>& gradient,
                     std::vector<double>& laplacian)
{
    std::vector<std::complex<double>> local_g(basis.npw);
    const FilledCavity& filled = response.solvent_aware;
    basis.real2recip(filled.local.data(), local_g.data());
    spec_derivs(local_g, basis, tpiba, gradient, laplacian);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        local_g[ig] *= probe_kernel[ig];
    }
    std::vector<ModuleBase::Vector3<double>> fraction_gradient;
    std::vector<double> fraction_laplacian;
    spec_derivs(local_g, basis, tpiba, fraction_gradient, fraction_laplacian);
    sa_chain_derivs(filled.local, filled.filling, fraction_gradient, fraction_laplacian, gradient,
                    laplacian);
}

// Environ eps_of_boundary: eps = exp(L (1 - s)), L = ln(eps_bulk), and
// the sqrt-CG coefficient eps (lapl ln eps / 2 + |grad ln eps|^2 / 4) / (4 pi).
void eps_of_boundary(const CavityParameters& cavity,
                     const std::vector<ModuleBase::Vector3<double>>& gradient,
                     const std::vector<double>& laplacian,
                     SccsResponse& response,
                     std::vector<double>& coefficient)
{
    const std::size_t size = response.solute.size();
    const double log_bulk = std::log(cavity.epsilon_bulk);
    coefficient.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        const double log_epsilon = log_bulk * (1.0 - response.solute[i]);
        response.epsilon[i] = std::exp(log_epsilon);
        const ModuleBase::Vector3<double>& local_gradient = gradient[i];
        const double gradient_square = local_gradient * local_gradient;
        response.grad_log_epsilon[i] = local_gradient * (-log_bulk);
        const double lap_log = -log_bulk * laplacian[i];
        coefficient[i] = response.epsilon[i] * (0.5 * lap_log + 0.25 * log_bulk * log_bulk * gradient_square)
                         / ModuleBase::FOUR_PI;
    }
}

// Environ solvent_aware_boundary: keep s as the local boundary and make solute
// the filled s_sa. The chain grad n, grad c = p * s' grad n and s'' serve the
// chain surface with every boundary condition; density_g is the transform of
// the cavity density.
void fill_cavity(const std::vector<double>& density,
                 const CavityParameters& cavity,
                 const std::vector<double>& probe_kernel,
                 const ModulePW::PW_Basis& basis,
                 double tpiba,
                 SccsResponse& response,
                 std::vector<std::complex<double>>& density_g)
{
    const std::size_t size = density.size();
    FilledCavity& filled = response.solvent_aware;
    filled.active = true;
    // Move s out and s_sa in without copies; filling.boundary is left empty.
    filled.local.swap(response.solute);
    filled.filling = sa_boundary(filled.local, probe_kernel, cavity.solvent_aware, basis);
    response.solute.swap(filled.filling.boundary);
    const std::vector<double>& boundary = response.solute;
    const std::vector<double>& local = filled.local;
    const auto add_filling = [&boundary, &local](std::size_t i, std::array<double, 1>& sum) {
        sum[0] += boundary[i] - local[i];
    };
    double filled_volume = thread_sums<1>(size, add_filling)[0];
    Parallel_Reduce::reduce_pool(filled_volume);
    filled.volume = filled_volume * basis.omega / basis.nxyz;

    density_g.resize(basis.npw);
    basis.real2recip(density.data(), density_g.data());
    filled.density_gradient.resize(size);
    XC_Functional::grad_rho(density_g.data(), filled.density_gradient.data(), &basis, tpiba);
    std::vector<ModuleBase::Vector3<double>> local_gradient(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        local_gradient[i] = filled.density_gradient[i] * response.dsolute_drho[i];
    }
    filled.fraction_gradient = convolve_grad(probe_kernel, basis, local_gradient);
}

// Environ deriv_method 'chain': grad s = s' grad n and
// lapl s = s' lapl n + s'' |grad n|^2 from the FFT derivatives of n.
void chain_derivs(const std::vector<std::complex<double>>& density_g,
                  const std::vector<ModuleBase::Vector3<double>>& density_gradient,
                  const std::vector<double>& dsolute_drho,
                  const std::vector<double>& d2solute_drho2,
                  const ModulePW::PW_Basis& basis,
                  double tpiba,
                  std::vector<ModuleBase::Vector3<double>>& gradient,
                  std::vector<double>& laplacian)
{
    const std::size_t size = density_gradient.size();
    std::vector<double> density_laplacian(size);
    // These XC helpers use G in units of tpiba and ABACUS's normalized FFT.
    XC_Functional::laplacian_rho(density_g.data(), density_laplacian.data(), &basis, tpiba);
    gradient.resize(size);
    laplacian.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        const double gradient_square = density_gradient[i] * density_gradient[i];
        gradient[i] = density_gradient[i] * dsolute_drho[i];
        laplacian[i] = dsolute_drho[i] * density_laplacian[i] + d2solute_drho2[i] * gradient_square;
    }
}
} // namespace

void make_cav_derivs(const std::vector<double>& density,
                     const CavityParameters& cavity,
                     const std::vector<double>& probe_kernel,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     bool open_boundary,
                     SccsResponse& response,
                     CavityDerivatives& derivatives)
{
    derivatives.filter = switch_filter(cavity, basis);
    const bool lowpass = uses_lowpass(cavity);
    // The configuration, not the local array sizes, decides the filling.
    const bool filled = uses_sa(cavity.solvent_aware);
    const std::size_t kernel_size = basis.npw;
    if (filled && probe_kernel.size() != kernel_size)
    {
        ModuleBase::WARNING_QUIT("ModuleSccs::make_cav_derivs",
                                 "solvent-aware SCCS needs the probe kernel on the local G vectors");
    }
    derivatives.gradient.clear();
    std::vector<double>& coefficient = derivatives.coefficient;
    const std::size_t size = density.size();
    response.solute.resize(size);
    response.dsolute_drho.resize(size);
    response.epsilon.resize(size);
    response.grad_log_epsilon.resize(size);
    response.solvent_aware = FilledCavity();
    std::vector<double> d2solute_drho2(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        const CavityPoint point = evaluate_cavity(density[i], cavity);
        response.solute[i] = point.solute;
        response.dsolute_drho[i] = point.dsolute_drho;
        response.epsilon[i] = point.epsilon;
        d2solute_drho2[i] = point.d2solute_drho2;
    }
    std::vector<std::complex<double>> density_g;
    if (filled)
    {
        response.solvent_aware.d2solute_drho2.swap(d2solute_drho2);
        fill_cavity(density, cavity, probe_kernel, basis, tpiba, response, density_g);
    }

    if (open_boundary)
    {
        // Open boundaries differentiate the switching boundary on the FFT grid,
        // since its chain coefficient is not converged at the cavity edge.
        // Unfiltered FFT derivatives of s_sa ring where the filling switches
        // within a grid spacing and stall the SCF, so without lowpass the
        // filling is applied analytically to those of the local s. With
        // lowpass, s_sa is differentiated on the filtered grid (Environ fft
        // with deriv_lowpass): the analytic filling would put unfiltered
        // pointwise f'' and f''' into the exact cavity derivative.
        std::vector<ModuleBase::Vector3<double>> gradient;
        std::vector<double> laplacian;
        if (filled && !lowpass)
        {
            fft_fill_derivs(probe_kernel, basis, tpiba, response, gradient, laplacian);
        }
        else
        {
            fft_derivatives(response.solute, derivatives.filter, basis, tpiba, gradient, laplacian);
        }
        eps_of_boundary(cavity, gradient, laplacian, response, coefficient);
        if (lowpass)
        {
            derivatives.gradient.swap(gradient);
        }
        return;
    }

    // The filled cavity already holds the transform and gradient of n and s''.
    std::vector<ModuleBase::Vector3<double>> plain_gradient;
    if (!filled)
    {
        density_g.resize(basis.npw);
        basis.real2recip(density.data(), density_g.data());
        plain_gradient.resize(size);
        XC_Functional::grad_rho(density_g.data(), plain_gradient.data(), &basis, tpiba);
    }
    const FilledCavity& filling = response.solvent_aware;
    const std::vector<ModuleBase::Vector3<double>>& density_gradient
        = filled ? filling.density_gradient : plain_gradient;
    const std::vector<double>& second_derivative = filled ? filling.d2solute_drho2 : d2solute_drho2;
    std::vector<ModuleBase::Vector3<double>> gradient;
    std::vector<double> laplacian;
    chain_derivs(density_g, density_gradient, response.dsolute_drho, second_derivative, basis, tpiba, gradient,
                 laplacian);
    if (filled)
    {
        // Environ deriv_method 'chain' through the filling.
        const std::vector<double> fraction_laplacian = convolve_probe(probe_kernel, basis, laplacian);
        sa_chain_derivs(filling.local, filling.filling, filling.fraction_gradient,
                        fraction_laplacian, gradient, laplacian);
    }
    eps_of_boundary(cavity, gradient, laplacian, response, coefficient);
}

} // namespace ModuleSccs
