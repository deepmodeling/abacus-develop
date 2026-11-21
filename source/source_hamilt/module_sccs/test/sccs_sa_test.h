#ifndef SCCS_SA_TEST_H
#define SCCS_SA_TEST_H

#include "sccs_test.h"
#include "../sccs_cavity.h"
#include "../sccs_sa.h"

#include <vector>

// Grid helpers and the cubic-cell fixture of the solvent-aware tests
// (definitions in sccs_sa_test.cpp).
namespace SccsTest
{
const ModuleBase::Matrix3 cubic_lattice(1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0);
const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);

void setup_basis(ModulePW::PW_Basis& basis, double scale, const ModuleBase::Matrix3& lattice, double ecut);

// Fold a displacement into the cube [-length/2, length/2).
void fold_into_cube(double length, ModuleBase::Vector3<double>& r);

// Minimum-image displacement of grid point ir from the origin in a cube.
ModuleBase::Vector3<double> cubic_displacement(const ModulePW::PW_Basis& basis, double length, int ir);

std::vector<ModuleBase::Vector3<double>> spectral_gradient(const std::vector<double>& values,
                                                           const ModulePW::PW_Basis& basis,
                                                           double tpiba);

std::vector<double> spectral_laplacian(const std::vector<double>& values,
                                       const ModulePW::PW_Basis& basis,
                                       double tpiba);

double pool_dot(const std::vector<double>& left, const std::vector<double>& right);

double pool_max(const ModulePW::PW_Basis& basis, double value);

// values + step * direction, for central differences.
std::vector<double> shifted(const std::vector<double>& values,
                            const std::vector<double>& direction,
                            double step);

ModuleSccs::SolventAwareParameters probe_parameters(double solvent_radius);

struct DensitySurface
{
    double surface = 0.0;
    ModuleSccs::SolventAwareNonelectrostatic chain;
};

// A cubic cell and the probe kernel of the solvent-aware parameters.
class SolventAwareTest : public testing::Test
{
protected:
    SolventAwareTest() : basis("cpu", "double") {}

    void set_up_cube(double edge, double ecut, const ModuleSccs::SolventAwareParameters& probe);

    ModuleBase::Vector3<double> displacement(int ir) const { return cubic_displacement(basis, length, ir); }

    // Smooth sphere erfc((|r - center| - radius)/spread)/2.
    std::vector<double> sphere(const ModuleBase::Vector3<double>& center, double radius, double spread) const;

    ModuleSccs::SolventAwareBoundary fill(const std::vector<double>& local) const
    {
        return ModuleSccs::sa_boundary(local, kernel, parameters, basis);
    }

    // Two Gaussians of density along x, centred at +-offset; their midpoint is a density saddle.
    std::vector<double> gaussian_pair(double offset) const;

    DensitySurface surface_of_density(const std::vector<double>& density,
                                      const ModuleSccs::CavityParameters& cavity,
                                      double regularization) const;

    // Central difference of the discrete surface along direction against the
    // projected analytic derivative; returns the relative error.
    double surface_derivative_error(const std::vector<double>& density,
                                    const std::vector<double>& direction,
                                    const ModuleSccs::CavityParameters& cavity,
                                    double regularization,
                                    double step) const;

    double length = 0.0;
    double tpiba = 0.0;
    ModulePW::PW_Basis basis;
    ModuleSccs::SolventAwareParameters parameters;
    std::vector<double> kernel;
};
} // namespace SccsTest

#endif
