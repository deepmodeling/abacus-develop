#include "sccs_lowpass.h"
#include "sccs_cavity.h"
#include "sccs_cav_deriv.h"
#include "sccs_response.h"

#include "source_base/constants.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <utility>

namespace ModuleSccs
{
std::vector<double> switch_filter(const CavityParameters& cavity, const ModulePW::PW_Basis& basis)
{
    const bool lowpass = uses_lowpass(cavity);
    std::vector<double> filter;
    if (lowpass)
    {
        filter.resize(basis.npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            // Environ deriv_lowpass: Gcut is the density cutoff, not ecutwfc.
            const double argument = cavity.lowpass_p1 * basis.gg[ig] / basis.ggecut - cavity.lowpass_p2;
            filter[ig] = 0.5 * std::erfc(argument);
        }
    }
    return filter;
}

void lowpass_bnd_pot(const std::vector<double>& charge,
                     const CavityParameters& cavity,
                     const CavityDerivatives& derivatives,
                     const ModulePW::PW_Basis& basis,
                     double tpiba,
                     SccsResponse& response)
{
    const std::size_t size = charge.size();
    const double log_bulk = std::log(cavity.epsilon_bulk);
    const std::vector<double>& potential = response.polarization.potential;
    std::vector<double> weight(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        weight[i] = response.epsilon[i] * potential[i] * potential[i] / (8.0 * ModuleBase::PI);
    }
    // A = sqrt(eps) C^-1 sqrt(eps) + F is symmetric. For E=q^T A^-1 q/2,
    // dE/ds = L/2 (q v + lapl b + L div(b grad s)), b=eps v^2/(8 pi).
    // The filtered derivative transposes share one inverse FFT.
    std::vector<std::complex<double>> transpose_g(basis.npw);
    basis.real2recip(weight.data(), transpose_g.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        transpose_g[ig] *= -tpiba * tpiba * basis.gg[ig] * derivatives.filter[ig];
    }
    std::vector<double> component(size);
    std::vector<std::complex<double>> component_g(basis.npw);
    for (int d = 0; d < 3; ++d)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            component[i] = weight[i] * derivatives.gradient[i][d];
        }
        basis.real2recip(component.data(), component_g.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            transpose_g[ig] += log_bulk * ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][d]
                               * component_g[ig] * derivatives.filter[ig];
        }
    }
    std::vector<double> transpose(size);
    basis.recip2real(transpose_g.data(), transpose.data());
    std::vector<double> candidate(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        candidate[i] = 0.5 * log_bulk * (charge[i] * potential[i] + transpose[i]);
    }
    response.boundary_potential = std::move(candidate);
}
} // namespace ModuleSccs
