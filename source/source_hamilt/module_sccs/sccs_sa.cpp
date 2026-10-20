#include "sccs_sa.h"
#include "sccs_cavity.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>

namespace ModuleSccs
{
namespace
{
// Environ tools_math erfc threshold: dsfunct2 and d2sfunct2 vanish beyond it.
const double filling_argument_cutoff = 6.0;
// Environ function_erfc cuts the probe at width + 5 spread.
const double probe_cutoff_spreads = 5.0;

// Fractional grid coordinate folded into [-1/2, 1/2).
double folded_fraction(int index, int count)
{
    double fraction = static_cast<double>(index) / static_cast<double>(count);
    if (fraction >= 0.5)
    {
        fraction -= 1.0;
    }
    return fraction;
}

// Periodic images n with |r0 + n a_k| <= cutoff for r0 in the folded cell:
// the cell height along a_k is the volume over the area of the other two.
int image_count(double cutoff,
                const ModuleBase::Vector3<double>& first,
                const ModuleBase::Vector3<double>& second,
                double volume)
{
    const double reach = cutoff * (first ^ second).norm() / volume + 0.5;
    return static_cast<int>(std::ceil(reach));
}
} // namespace

std::vector<double> sa_probe_kernel(const ModulePW::PW_Basis& basis,
                                    const ModuleBase::Matrix3& lattice,
                                    double lattice_scale,
                                    const SolventAwareParameters& parameters)
{
    ModuleBase::timer::start("ModuleSccs", "sa_probe_kernel");
    const ModuleBase::Vector3<double> row1(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> row2(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> row3(lattice.e31, lattice.e32, lattice.e33);
    const ModuleBase::Vector3<double> a1 = row1 * lattice_scale;
    const ModuleBase::Vector3<double> a2 = row2 * lattice_scale;
    const ModuleBase::Vector3<double> a3 = row3 * lattice_scale;
    const double triple_product = a1 * (a2 ^ a3);
    const double volume = std::abs(triple_product);
    const double width = parameters.solvent_radius * parameters.radial_scale;
    const double spread = parameters.radial_spread;
    const double cutoff = width + probe_cutoff_spreads * spread;
    const int count_1 = image_count(cutoff, a2, a3, volume);
    const int count_2 = image_count(cutoff, a3, a1, volume);
    const int count_3 = image_count(cutoff, a1, a2, volume);

    std::vector<double> probe(basis.nrxx, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        const ModuleBase::Vector3<double> folded = a1 * folded_fraction(ix, basis.nx)
                                                   + a2 * folded_fraction(iy, basis.ny)
                                                   + a3 * folded_fraction(iz, basis.nz);
        double value = 0.0;
        for (int n1 = -count_1; n1 <= count_1; ++n1)
        {
            for (int n2 = -count_2; n2 <= count_2; ++n2)
            {
                for (int n3 = -count_3; n3 <= count_3; ++n3)
                {
                    const ModuleBase::Vector3<double> image = folded + a1 * static_cast<double>(n1)
                                                              + a2 * static_cast<double>(n2)
                                                              + a3 * static_cast<double>(n3);
                    const double distance = image.norm();
                    if (distance <= cutoff)
                    {
                        const double argument = (distance - width) / spread;
                        value += std::erfc(argument);
                    }
                }
            }
        }
        probe[ir] = value;
    }
    // The image sums dominate; a serial sum keeps the normalization independent
    // of the thread count.
    double integral = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        integral += probe[ir];
    }
    Parallel_Reduce::reduce_pool(integral);
    integral *= volume / static_cast<double>(basis.nxyz);
    // Environ convolution: fwfft(f1) fwfft(f2) omega with the 1/N forward
    // transform, so the kernel is omega times the transform of the normalized
    // probe. The probe is even on the grid, so its transform is real.
    std::vector<std::complex<double>> probe_g(basis.npw);
    basis.real2recip(probe.data(), probe_g.data());
    std::vector<double> kernel(basis.npw);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        kernel[ig] = probe_g[ig].real() * volume / integral;
    }
    ModuleBase::timer::end("ModuleSccs", "sa_probe_kernel");
    return kernel;
}

std::vector<double> convolve_probe(const std::vector<double>& kernel,
                                   const ModulePW::PW_Basis& basis,
                                   const std::vector<double>& values)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        values_g[ig] *= kernel[ig];
    }
    std::vector<double> result(values.size());
    basis.recip2real(values_g.data(), result.data());
    return result;
}

std::vector<ModuleBase::Vector3<double>> convolve_grad(
    const std::vector<double>& kernel,
    const ModulePW::PW_Basis& basis,
    const std::vector<ModuleBase::Vector3<double>>& gradient)
{
    const std::size_t size = gradient.size();
    std::vector<double> component(size);
    std::vector<ModuleBase::Vector3<double>> result(size);
    for (int d = 0; d < 3; ++d)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            component[i] = gradient[i][d];
        }
        const std::vector<double> convolved = convolve_probe(kernel, basis, component);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            result[i][d] = convolved[i];
        }
    }
    return result;
}

SolventAwareBoundary sa_boundary(const std::vector<double>& local,
                                 const std::vector<double>& kernel,
                                 const SolventAwareParameters& parameters,
                                 const ModulePW::PW_Basis& basis)
{
    SolventAwareBoundary filled;
    filled.fraction = convolve_probe(kernel, basis, local);
    const std::size_t size = local.size();
    const double threshold = parameters.filling_threshold;
    const double spread = parameters.filling_spread;
    const double sqrt_pi = std::sqrt(ModuleBase::PI);
    filled.filling.resize(size);
    filled.dfilling.assign(size, 0.0);
    filled.d2filling.assign(size, 0.0);
    filled.boundary.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        // Environ: fill = 1 - sfunct2, dfill = -dsfunct2, d2fill = -d2sfunct2.
        const double argument = (filled.fraction[i] - threshold) / spread;
        filled.filling[i] = 1.0 - 0.5 * std::erfc(argument);
        if (std::abs(argument) <= filling_argument_cutoff)
        {
            const double exponent = -argument * argument;
            const double gaussian = std::exp(exponent) / sqrt_pi;
            filled.dfilling[i] = gaussian / spread;
            filled.d2filling[i] = -2.0 * argument * gaussian / (spread * spread);
        }
        filled.boundary[i] = local[i] + (1.0 - local[i]) * filled.filling[i];
    }
    return filled;
}

void sa_chain_derivs(const std::vector<double>& local,
                     const SolventAwareBoundary& filled,
                     const std::vector<ModuleBase::Vector3<double>>& fraction_gradient,
                     const std::vector<double>& fraction_laplacian,
                     std::vector<ModuleBase::Vector3<double>>& gradient,
                     std::vector<double>& laplacian)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < local.size(); ++i)
    {
        const double empty = 1.0 - filled.filling[i];
        const double solvent = 1.0 - local[i];
        const double projection = gradient[i] * fraction_gradient[i];
        laplacian[i] = laplacian[i] * empty - 2.0 * projection * filled.dfilling[i]
                       + solvent * (filled.d2filling[i] * fraction_gradient[i].norm2()
                                    + filled.dfilling[i] * fraction_laplacian[i]);
        gradient[i] = gradient[i] * empty + fraction_gradient[i] * (solvent * filled.dfilling[i]);
    }
}

SolventAwareNonelectrostatic sa_nonelec(const FilledCavity& cavity,
                                        const std::vector<double>& dsolute_drho,
                                        const std::vector<double>& kernel,
                                        const ModulePW::PW_Basis& basis,
                                        double tpiba,
                                        double regularization,
                                        double surface_tension,
                                        double pressure)
{
    ModuleBase::timer::start("ModuleSccs", "sa_nonelec");
    const std::vector<double>& local = cavity.local;
    const SolventAwareBoundary& filled = cavity.filling;
    const std::vector<ModuleBase::Vector3<double>>& density_gradient = cavity.density_gradient;
    const std::vector<double>& d2solute_drho2 = cavity.d2solute_drho2;
    const std::vector<ModuleBase::Vector3<double>>& fraction_gradient = cavity.fraction_gradient;
    const std::size_t size = local.size();
    const double regularization_square = regularization * regularization;
    SolventAwareNonelectrostatic result;
    result.gradient.resize(size);
    std::vector<ModuleBase::Vector3<double>> filled_part(size);
    std::vector<ModuleBase::Vector3<double>> fraction_part(size);
    std::vector<double> convolved(size);
    std::vector<double> pointwise(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        const ModuleBase::Vector3<double> local_gradient = density_gradient[i] * dsolute_drho[i];
        const double empty = 1.0 - filled.filling[i];
        const double solvent = 1.0 - local[i];
        const double fraction_weight = solvent * filled.dfilling[i];
        const ModuleBase::Vector3<double> g = local_gradient * empty + fraction_gradient[i] * fraction_weight;
        result.gradient[i] = g;
        const double norm_square = g.norm2() + regularization_square;
        const ModuleBase::Vector3<double> unit = g / std::sqrt(norm_square);
        const double local_projection = unit * local_gradient;
        const double fraction_projection = unit * fraction_gradient[i];
        filled_part[i] = unit * empty;
        fraction_part[i] = unit * fraction_weight;
        // gamma alpha + p (1 - s) f' and gamma beta + p (1 - f): the convolved
        // and pointwise coefficients of the variation of s.
        const double alpha = -filled.dfilling[i] * local_projection
                             + solvent * filled.d2filling[i] * fraction_projection;
        const double beta = -filled.dfilling[i] * fraction_projection;
        convolved[i] = surface_tension * alpha + pressure * fraction_weight;
        pointwise[i] = surface_tension * beta + pressure * empty;
    }
    // A = (1 - f) u + p * ((1 - s) f' u): the coefficient of the variation of a.
    const std::vector<ModuleBase::Vector3<double>> fraction_adjoint = convolve_grad(kernel, basis,
                                                                                   fraction_part);
    const std::vector<double> convolved_adjoint = convolve_probe(kernel, basis, convolved);
    std::vector<double> component(size);
    std::vector<std::complex<double>> component_g(basis.npw);
    std::vector<std::complex<double>> divergence_g(basis.npw, std::complex<double>(0.0, 0.0));
    result.density_potential.resize(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        filled_part[i] += fraction_adjoint[i];
        const double curvature = d2solute_drho2[i] * (filled_part[i] * density_gradient[i]);
        result.density_potential[i] = dsolute_drho[i] * (pointwise[i] + convolved_adjoint[i])
                                      + surface_tension * curvature;
    }
    // a = s' D n with the spectral gradient D, whose transpose is -div.
    for (int d = 0; d < 3; ++d)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < size; ++i)
        {
            component[i] = dsolute_drho[i] * filled_part[i][d];
        }
        basis.real2recip(component.data(), component_g.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            divergence_g[ig] += ModuleBase::IMAG_UNIT * tpiba * basis.gcar[ig][d] * component_g[ig];
        }
    }
    std::vector<double> divergence(size);
    basis.recip2real(divergence_g.data(), divergence.data());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        result.density_potential[i] -= surface_tension * divergence[i];
    }
    ModuleBase::timer::end("ModuleSccs", "sa_nonelec");
    return result;
}

std::vector<double> sa_adjoint(const FilledCavity& cavity,
                               const std::vector<double>& kernel,
                               const ModulePW::PW_Basis& basis,
                               const std::vector<double>& boundary_potential)
{
    const std::vector<double>& local = cavity.local;
    const SolventAwareBoundary& filled = cavity.filling;
    const std::size_t size = local.size();
    std::vector<double> weighted(size);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        weighted[i] = (1.0 - local[i]) * filled.dfilling[i] * boundary_potential[i];
    }
    std::vector<double> result = convolve_probe(kernel, basis, weighted);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        result[i] += (1.0 - filled.filling[i]) * boundary_potential[i];
    }
    return result;
}
} // namespace ModuleSccs
