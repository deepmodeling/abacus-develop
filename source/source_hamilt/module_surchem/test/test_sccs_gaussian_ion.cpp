#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_gaussian_ion.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace
{

// One ion of valence zv at fractional position tau in a 10 bohr cubic cell.
void setup_single_ion(ModulePW::PW_Basis& basis,
                      UnitCell& cell,
                      const double zv,
                      const ModuleBase::Vector3<double>& tau)
{
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 80.0);
    basis.initparameters(false, 80.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = length * length * length;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].ncpp.zv = zv;
    cell.atoms[0].tau.push_back(tau);
}

TEST(SccsGaussianIon, UsesEnvironSpreadConvention)
{
    ModulePW::PW_Basis basis("cpu", "double");
    UnitCell cell;
    const ModuleBase::Vector3<double> tau(0.37, 0.43, 0.52);
    setup_single_ion(basis, cell, 2.0, tau);
    const double spread = ModuleSccs::gaussian_ion_spread;
    EXPECT_DOUBLE_EQ(spread, 0.5);

    const std::vector<double> density = ModuleSccs::gaussian_ionic_density(cell, basis, spread);
    ASSERT_EQ(density.size(), static_cast<std::size_t>(basis.nrxx));
    // ENVIRON atomicspread: rho(r) = zv exp(-r^2/s^2) / (sqrt(pi) s)^3, so that
    // rho(G) = zv/Omega exp(-s^2 G^2 / 4).
    std::vector<std::complex<double>> density_g(basis.npw);
    basis.real2recip(density.data(), density_g.data());
    double maximum_error = 0.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const double g2 = cell.tpiba2 * basis.gg[ig];
        const double expected = 2.0 / cell.omega * std::exp(-0.25 * spread * spread * g2);
        const double error = std::abs(std::abs(density_g[ig]) - expected);
        maximum_error = std::max(maximum_error, error);
    }
    EXPECT_LT(maximum_error, 1.0e-14);

    double charge = 0.0;
    for (const double value : density)
    {
        charge += value;
    }
    charge *= cell.omega / static_cast<double>(basis.nxyz);
    EXPECT_NEAR(charge, 2.0, 1.0e-12);
}

TEST(SccsGaussianIon, ForceMatchesOpenQuadraticPotentialDerivative)
{
    ModulePW::PW_Basis basis("cpu", "double");
    UnitCell cell;
    const ModuleBase::Vector3<double> tau(0.37, 0.43, 0.52);
    setup_single_ion(basis, cell, 2.0, tau);
    const double length = cell.lat0;
    const double spread = ModuleSccs::gaussian_ion_spread;
    const double volume_element = cell.omega / basis.nxyz;
    // A minimum-image quadratic, like the PCC part of an open-boundary potential.
    std::vector<double> shape_potential(basis.nrxx);
    for (int index = 0; index < basis.nrxx; ++index)
    {
        const int iy = (index / basis.nplane) % basis.ny;
        const double y = length * iy / basis.ny - 0.5 * length;
        shape_potential[index] = 0.3 * y * y;
    }
    const ModuleBase::matrix force
        = ModuleSccs::gaussian_ionic_force(cell, basis, spread, shape_potential);
    const double displacement = 1.0e-5;
    cell.atoms[0].tau[0].y += displacement / length;
    const std::vector<double> plus = ModuleSccs::gaussian_ionic_density(cell, basis, spread);
    cell.atoms[0].tau[0].y -= 2.0 * displacement / length;
    const std::vector<double> minus = ModuleSccs::gaussian_ionic_density(cell, basis, spread);
    double energy_difference = 0.0;
    for (int index = 0; index < basis.nrxx; ++index)
    {
        energy_difference += (plus[index] - minus[index]) * shape_potential[index] * volume_element;
    }
    const double finite_difference = -energy_difference / (2.0 * displacement);
    EXPECT_NEAR(force(0, 1), finite_difference, 1.0e-8);
    EXPECT_NEAR(force(0, 0), 0.0, 1.0e-12);
    EXPECT_NEAR(force(0, 2), 0.0, 1.0e-12);
}

TEST(SccsGaussianIon, RejectsMismatchedPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
    UnitCell cell;
    const ModuleBase::Vector3<double> tau(0.5, 0.5, 0.5);
    setup_single_ion(basis, cell, 1.0, tau);
    const std::vector<double> wrong_size(basis.nrxx + 1, 0.0);
    EXPECT_THROW(ModuleSccs::gaussian_ionic_force(cell, basis, ModuleSccs::gaussian_ion_spread,
                                                  wrong_size),
                 std::runtime_error);
}

// A sulfur-like atom (zv 6) and a hydrogen (zv 1) in a 10 bohr cubic cell.
void setup_sulfur_hydrogen(ModulePW::PW_Basis& basis, UnitCell& cell)
{
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 80.0);
    basis.initparameters(false, 80.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = length * length * length;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 2;
    cell.nat = 2;
    cell.atoms = new Atom[2];
    cell.atoms[0].na = 1;
    cell.atoms[0].ncpp.zv = 6.0;
    cell.atoms[0].ncpp.psd = "S";
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.37, 0.43, 0.52));
    cell.atoms[1].na = 1;
    cell.atoms[1].ncpp.zv = 1.0;
    cell.atoms[1].ncpp.psd = " H ";
    cell.atoms[1].tau.push_back(ModuleBase::Vector3<double>(0.51, 0.43, 0.52));
}

TEST(SccsGaussianIon, CoreDensitySkipsHydrogen)
{
    ModulePW::PW_Basis basis("cpu", "double");
    UnitCell cell;
    setup_sulfur_hydrogen(basis, cell);
    const std::vector<double> core = ModuleSccs::gaussian_core_density(cell, basis, 0.5);
    const std::vector<double> ionic
        = ModuleSccs::gaussian_ionic_density(cell, basis, ModuleSccs::gaussian_ion_spread);
    const double volume_element = cell.omega / static_cast<double>(basis.nxyz);
    double core_charge = 0.0;
    double ionic_charge = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        core_charge += core[ir] * volume_element;
        ionic_charge += ionic[ir] * volume_element;
    }
    EXPECT_NEAR(core_charge, 6.0, 1.0e-10);
    EXPECT_NEAR(ionic_charge, 7.0, 1.0e-10);
}

TEST(SccsGaussianIon, CoreForceMatchesDerivativeAndSkipsHydrogen)
{
    ModulePW::PW_Basis basis("cpu", "double");
    UnitCell cell;
    setup_sulfur_hydrogen(basis, cell);
    const double length = cell.lat0;
    const double spread = 0.7;
    const double volume_element = cell.omega / basis.nxyz;
    std::vector<double> potential(basis.nrxx);
    for (int index = 0; index < basis.nrxx; ++index)
    {
        const int iy = (index / basis.nplane) % basis.ny;
        const double y = length * iy / basis.ny - 0.5 * length;
        potential[index] = 0.3 * y * y;
    }
    const ModuleBase::matrix force = ModuleSccs::gaussian_core_force(cell, basis, spread, potential);
    const double displacement = 1.0e-5;
    cell.atoms[0].tau[0].y += displacement / length;
    const std::vector<double> plus = ModuleSccs::gaussian_core_density(cell, basis, spread);
    cell.atoms[0].tau[0].y -= 2.0 * displacement / length;
    const std::vector<double> minus = ModuleSccs::gaussian_core_density(cell, basis, spread);
    cell.atoms[0].tau[0].y += displacement / length;
    double energy_difference = 0.0;
    for (int index = 0; index < basis.nrxx; ++index)
    {
        energy_difference += (plus[index] - minus[index]) * potential[index] * volume_element;
    }
    const double finite_difference = -energy_difference / (2.0 * displacement);
    EXPECT_NEAR(force(0, 1), finite_difference, 1.0e-7);
    EXPECT_GT(std::abs(force(0, 1)), 1.0e-3);
    for (int d = 0; d < 3; ++d)
    {
        EXPECT_DOUBLE_EQ(force(1, d), 0.0);
    }
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
