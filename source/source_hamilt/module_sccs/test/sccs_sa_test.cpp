#include "sccs_sa_test.h"

#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>

namespace SccsTest
{
void setup_basis(ModulePW::PW_Basis& basis, double scale, const ModuleBase::Matrix3& lattice, double ecut)
{
#ifdef __MPI
    basis.initmpi(pool_size, pool_rank, POOL_WORLD);
#endif
    basis.initgrids(scale, lattice, ecut);
    basis.initparameters(false, ecut, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
}

void fold_into_cube(double length, ModuleBase::Vector3<double>& r)
{
    for (int d = 0; d < 3; ++d)
    {
        const double fraction = r[d] / length;
        r[d] -= length * std::round(fraction);
    }
}

ModuleBase::Vector3<double> cubic_displacement(const ModulePW::PW_Basis& basis, double length, int ir)
{
    const int ix = ir / (basis.ny * basis.nplane);
    const int iy = ir / basis.nplane - ix * basis.ny;
    const int iz = ir % basis.nplane + basis.startz_current;
    const double x = length * ix / basis.nx;
    const double y = length * iy / basis.ny;
    const double z = length * iz / basis.nz;
    ModuleBase::Vector3<double> r(x, y, z);
    fold_into_cube(length, r);
    return r;
}

std::vector<ModuleBase::Vector3<double>> spectral_gradient(const std::vector<double>& values,
                                                           const ModulePW::PW_Basis& basis,
                                                           double tpiba)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    std::vector<ModuleBase::Vector3<double>> gradient(values.size());
    XC_Functional::grad_rho(values_g.data(), gradient.data(), &basis, tpiba);
    return gradient;
}

std::vector<double> spectral_laplacian(const std::vector<double>& values,
                                       const ModulePW::PW_Basis& basis,
                                       double tpiba)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    std::vector<double> laplacian(values.size());
    XC_Functional::laplacian_rho(values_g.data(), laplacian.data(), &basis, tpiba);
    return laplacian;
}

double pool_dot(const std::vector<double>& left, const std::vector<double>& right)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        sum += left[i] * right[i];
    }
    Parallel_Reduce::reduce_pool(sum);
    return sum;
}

double pool_max(const ModulePW::PW_Basis& basis, double value)
{
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, value);
    return value;
}

std::vector<double> shifted(const std::vector<double>& values,
                            const std::vector<double>& direction,
                            double step)
{
    std::vector<double> result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        result[i] = values[i] + step * direction[i];
    }
    return result;
}

ModuleSccs::SolventAwareParameters probe_parameters(double solvent_radius)
{
    ModuleSccs::SolventAwareParameters parameters;
    parameters.solvent_radius = solvent_radius;
    return parameters;
}

void SolventAwareTest::set_up_cube(double edge, double ecut, const ModuleSccs::SolventAwareParameters& probe)
{
    length = edge;
    tpiba = ModuleBase::TWO_PI / edge;
    parameters = probe;
    setup_basis(basis, edge, cubic_lattice, ecut);
    kernel = ModuleSccs::sa_probe_kernel(basis, cubic_lattice, edge, probe);
}

std::vector<double> SolventAwareTest::sphere(const ModuleBase::Vector3<double>& center, double radius, double spread) const
{
    std::vector<double> values(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        ModuleBase::Vector3<double> r = displacement(ir) - center;
        fold_into_cube(length, r);
        const double argument = (r.norm() - radius) / spread;
        values[ir] = 0.5 * std::erfc(argument);
    }
    return values;
}

std::vector<double> SolventAwareTest::gaussian_pair(double offset) const
{
    std::vector<double> density(basis.nrxx);
    const ModuleBase::Vector3<double> first(offset, 0.0, 0.0);
    const ModuleBase::Vector3<double> second(-offset, 0.0, 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const ModuleBase::Vector3<double> r = displacement(ir);
        const double exponent_first = -(r - first).norm2();
        const double exponent_second = -(r - second).norm2();
        const double first_gaussian = std::exp(exponent_first);
        const double second_gaussian = std::exp(exponent_second);
        density[ir] = 4.8e-3 * (first_gaussian + second_gaussian);
    }
    return density;
}

DensitySurface SolventAwareTest::surface_of_density(const std::vector<double>& density,
                                                    const ModuleSccs::CavityParameters& cavity,
                                                    double regularization) const
{
    const std::size_t size = density.size();
    ModuleSccs::FilledCavity filled;
    filled.active = true;
    filled.local.resize(size);
    filled.d2solute_drho2.resize(size);
    std::vector<double> dsolute(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(density[i], cavity);
        filled.local[i] = point.solute;
        dsolute[i] = point.dsolute_drho;
        filled.d2solute_drho2[i] = point.d2solute_drho2;
    }
    filled.filling = fill(filled.local);
    filled.density_gradient = spectral_gradient(density, basis, tpiba);
    std::vector<ModuleBase::Vector3<double>> local_gradient(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        local_gradient[i] = filled.density_gradient[i] * dsolute[i];
    }
    filled.fraction_gradient = ModuleSccs::convolve_grad(kernel, basis, local_gradient);
    // Unit surface tension and no pressure: density_potential is dS/dn.
    const double surface_tension = 1.0;
    const double pressure = 0.0;
    DensitySurface result;
    result.chain = ModuleSccs::sa_nonelec(filled, dsolute, kernel, basis, tpiba, regularization,
                                          surface_tension, pressure);
    for (std::size_t i = 0; i < size; ++i)
    {
        const double norm_square = result.chain.gradient[i].norm2() + regularization * regularization;
        result.surface += std::sqrt(norm_square) - regularization;
    }
    Parallel_Reduce::reduce_pool(result.surface);
    return result;
}

double SolventAwareTest::surface_derivative_error(const std::vector<double>& density,
                                                  const std::vector<double>& direction,
                                                  const ModuleSccs::CavityParameters& cavity,
                                                  double regularization,
                                                  double step) const
{
    const DensitySurface center = surface_of_density(density, cavity, regularization);
    const double backward = -step;
    const std::vector<double> plus = shifted(density, direction, step);
    const std::vector<double> minus = shifted(density, direction, backward);
    const double surface_plus = surface_of_density(plus, cavity, regularization).surface;
    const double surface_minus = surface_of_density(minus, cavity, regularization).surface;
    const double finite_difference = (surface_plus - surface_minus) / (2.0 * step);
    const double analytic = pool_dot(center.chain.density_potential, direction);
    EXPECT_GT(std::abs(finite_difference), 1.0e-3);
    const double miss = std::abs(analytic - finite_difference);
    return miss / std::abs(finite_difference);
}

} // namespace SccsTest
