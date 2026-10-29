#include "sccs_cav_deriv.h"
#include "sccs_cavity.h"
#include "sccs_lowpass.h"
#include "sccs_response.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
namespace
{
// Environ eps_of_boundary: eps = exp(L (1 - s)), L = ln(eps_bulk), and
// the sqrt-CG coefficient eps (lapl ln eps / 2 + |grad ln eps|^2 / 4) / (4 pi)
// from grad s and lapl s.
void eps_of_boundary(const CavityParameters& cavity,
                     const std::vector<ModuleBase::Vector3<double>>& gradient,
                     const std::vector<double>& laplacian,
                     SccsResponse& response,
                     std::vector<double>& coefficient)
{
    const std::size_t size = response.solute.size();
    const double log_bulk = std::log(cavity.epsilon_bulk);
    coefficient.resize(size);
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
} // namespace

void make_cav_derivs(const std::vector<double>& density,
                     const CavityParameters& cavity,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     bool open_boundary,
                     SccsResponse& response,
                     CavityDerivatives& derivatives)
{
    derivatives.filter = switch_filter(cavity, basis);
    const bool lowpass = uses_lowpass(cavity);
    derivatives.gradient.clear();
    std::vector<double>& coefficient = derivatives.coefficient;
    const std::size_t size = density.size();
    response.solute.resize(size);
    response.dsolute_drho.resize(size);
    response.epsilon.resize(size);
    response.depsilon_drho.resize(size);
    response.grad_log_epsilon.resize(size);
    std::vector<double> d2solute_drho2(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const CavityPoint point = evaluate_cavity(density[i], cavity);
        response.solute[i] = point.solute;
        response.dsolute_drho[i] = point.dsolute_drho;
        response.epsilon[i] = point.epsilon;
        response.depsilon_drho[i] = point.depsilon_drho;
        d2solute_drho2[i] = point.d2solute_drho2;
    }

    std::vector<ModuleBase::Vector3<double>> gradient(size);
    std::vector<double> laplacian(size);
    if (open_boundary)
    {
        // Original PCC path: differentiate the switching boundary on the FFT
        // grid, since its chain coefficient is not converged at the cavity edge.
        std::vector<std::complex<double>> boundary_g(basis.npw);
        basis.real2recip(response.solute.data(), boundary_g.data());
        if (lowpass)
        {
            for (int ig = 0; ig < basis.npw; ++ig)
            {
                boundary_g[ig] *= derivatives.filter[ig];
            }
        }
        XC_Functional::grad_rho(boundary_g.data(), gradient.data(), &basis, tpiba);
        XC_Functional::laplacian_rho(boundary_g.data(), laplacian.data(), &basis, tpiba);
    }
    else
    {
        // Environ deriv_method 'chain': grad s = s' grad n and
        // lapl s = s' lapl n + s'' |grad n|^2 from the FFT derivatives of n.
        std::vector<std::complex<double>> density_g(basis.npw);
        basis.real2recip(density.data(), density_g.data());
        std::vector<ModuleBase::Vector3<double>> density_gradient(size);
        std::vector<double> density_laplacian(size);
        // These XC helpers use G in units of tpiba and ABACUS's normalized FFT.
        XC_Functional::grad_rho(density_g.data(), density_gradient.data(), &basis, tpiba);
        XC_Functional::laplacian_rho(density_g.data(), density_laplacian.data(), &basis, tpiba);
        for (std::size_t i = 0; i < size; ++i)
        {
            const double gradient_square = density_gradient[i] * density_gradient[i];
            gradient[i] = density_gradient[i] * response.dsolute_drho[i];
            laplacian[i] = response.dsolute_drho[i] * density_laplacian[i] + d2solute_drho2[i] * gradient_square;
        }
    }
    eps_of_boundary(cavity, gradient, laplacian, response, coefficient);
    if (open_boundary && lowpass)
    {
        derivatives.gradient.swap(gradient);
    }
}

} // namespace ModuleSccs
