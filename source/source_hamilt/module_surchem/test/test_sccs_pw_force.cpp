#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_pw_charge.h"
#include "../sccs/sccs_pw_force.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

namespace
{

double ionic_coupling_energy(const UnitCell& cell,
                             const ModulePW::PW_Basis& basis,
                             const ModuleBase::matrix& radial_local_potential,
                             const std::vector<double>& reaction_potential)
{
    std::vector<std::complex<double>> local_potential_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const std::complex<double> phase
            = std::exp(ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI
                       * (basis.gcar[ig] * cell.atoms[0].tau[0]));
        local_potential_g[ig] = radial_local_potential(0, basis.ig2igg[ig]) * phase;
    }
    std::vector<double> local_potential(basis.nrxx);
    basis.recip2real(local_potential_g.data(), local_potential.data());
    const std::vector<double> ionic_density
        = ModuleSccs::ionic_charge_from_local_potential(local_potential,
                                                        cell.atoms[0].ncpp.zv,
                                                        cell.omega,
                                                        cell.tpiba,
                                                        basis);
    double energy = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        energy += ionic_density[ir] * reaction_potential[ir];
    }
    return energy * cell.omega / static_cast<double>(basis.nxyz);
}

TEST(SccsPwForce, MatchesTranslatedSmoothChargeFiniteDifference)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 30.0);
    basis.initparameters(false, 30.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.collect_uniqgg();

    UnitCell cell;
    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = length * length * length;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.37, 0.43, 0.52));

    ModuleBase::matrix radial_local_potential(1, basis.ngg);
    const double gaussian_width = 0.35;
    for (int radial_index = 0; radial_index < basis.ngg; ++radial_index)
    {
        const double gg = basis.gg_uniq[radial_index];
        if (gg == 0.0)
        {
            radial_local_potential(0, radial_index) = 0.0;
            continue;
        }
        const double coulomb_rydberg
            = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * gg);
        radial_local_potential(0, radial_index)
            = -coulomb_rydberg / cell.omega
              * std::exp(-gaussian_width * cell.tpiba2 * gg);
    }

    std::vector<double> reaction_potential(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const double fractional_x
            = (static_cast<double>(ix) + 0.5) / static_cast<double>(basis.nx);
        reaction_potential[ir]
            = 0.7 * std::sin(ModuleBase::TWO_PI * fractional_x)
              + 0.2 * std::cos(2.0 * ModuleBase::TWO_PI * fractional_x);
    }

    const ModuleBase::matrix force
        = ModuleSccs::smooth_ionic_force_hartree(cell,
                                                 basis,
                                                 radial_local_potential,
                                                 reaction_potential);
    const double displacement = 1.0e-4;
    cell.atoms[0].tau[0].x += displacement / length;
    const double energy_plus
        = ionic_coupling_energy(cell, basis, radial_local_potential, reaction_potential);
    cell.atoms[0].tau[0].x -= 2.0 * displacement / length;
    const double energy_minus
        = ionic_coupling_energy(cell, basis, radial_local_potential, reaction_potential);
    cell.atoms[0].tau[0].x += displacement / length;
    const double finite_difference_force
        = -(energy_plus - energy_minus) / (2.0 * displacement);

    EXPECT_NEAR(force(0, 0), finite_difference_force, 1.0e-9);
    EXPECT_NEAR(force(0, 1), 0.0, 1.0e-12);
    EXPECT_NEAR(force(0, 2), 0.0, 1.0e-12);
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
