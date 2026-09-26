#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../surchem.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace
{

std::vector<double> local_potential_for_atom(const ModulePW::PW_Basis& basis,
                                             const ModuleBase::matrix& radial_local_potential,
                                             const ModuleBase::Vector3<double>& tau)
{
    std::vector<std::complex<double>> potential_g(basis.npw);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const std::complex<double> phase
            = std::exp(ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI
                       * (basis.gcar[ig] * tau));
        potential_g[ig] = radial_local_potential(0, basis.ig2igg[ig]) * phase;
    }
    std::vector<double> potential(basis.nrxx);
    basis.recip2real(potential_g.data(), potential.data());
    return potential;
}

double evaluate_reaction_energy(surchem& solvent,
                                const UnitCell& cell,
                                const ModulePW::PW_Basis& basis,
                                const ModuleBase::matrix& radial_local_potential,
                                const std::vector<double>& electron_density)
{
    const std::vector<double> local_potential
        = local_potential_for_atom(basis,
                                   radial_local_potential,
                                   cell.atoms[0].tau[0]);
    const double* density_channels[1] = {electron_density.data()};
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);
    return solvent.sccs_result().electrostatic.reaction_energy;
}

double evaluate_total_electrostatic_energy(
    surchem& solvent,
    const UnitCell& cell,
    const ModulePW::PW_Basis& basis,
    const ModuleBase::matrix& radial_local_potential,
    const std::vector<double>& electron_density)
{
    evaluate_reaction_energy(solvent,
                             cell,
                             basis,
                             radial_local_potential,
                             electron_density);
    return solvent.sccs_result().electrostatic.reaction_energy
           + solvent.sccs_result().vacuum_pcc_energy;
}

TEST(SolForce, ConvertsPointIonPccForceFromHartreeToRydberg)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    const double volume = length * length * length;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.collect_uniqgg();

    UnitCell cell;
    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = volume;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 2;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 2;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 0.5;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.4, 0.5, 0.5));
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.6, 0.5, 0.5));

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 0.0;
    parameters.expected_ionic_charge = 1.0;
    parameters.sccs_config.cavity.density_min = 1.0e-2;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 1.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc0d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc0d;
    parameters.sccs_config.max_iterations = 20;
    parameters.sccs_config.mixing = 0.7;
    parameters.sccs_config.tolerance_rms = 1.0e-14;
    parameters.sccs_config.tolerance_max = 1.0e-14;

    surchem solvent;
    solvent.set_parameters(parameters);
    std::vector<double> electron_density(basis.nrxx, 0.0);
    const double* density_channels[1] = {electron_density.data()};
    std::vector<double> local_potential(basis.nrxx, 0.0);
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);

    ModuleBase::matrix unused_vloc(1, basis.ngg);
    ModuleBase::matrix force(2, 3);
    solvent.cal_force_sol(cell, &basis, unused_vloc, 1, force);

    ModuleSccs::PccGeometry geometry
        = ModuleSccs::pcc_geometry(lattice, length, 1.0e-10);
    for (int atom = 0; atom < 2; ++atom)
    {
        ModuleSccs::PointCharge point;
        point.charge = 0.5;
        point.position = cell.atoms[0].tau[atom] * length;
        const ModuleBase::Vector3<double> expected_hartree
            = ModuleSccs::pcc_point_charge_force(solvent.sccs_result().point_solute_moments,
                                                 point,
                                                 geometry);
        EXPECT_NEAR(force(atom, 0), 2.0 * expected_hartree.x, 1.0e-14);
        EXPECT_NEAR(force(atom, 1), 2.0 * expected_hartree.y, 1.0e-14);
        EXPECT_NEAR(force(atom, 2), 2.0 * expected_hartree.z, 1.0e-14);
    }
    EXPECT_NEAR(force(0, 0) + force(1, 0), 0.0, 1.0e-14);
}

TEST(SolForce, RequiresPwBasis)
{
    UnitCell cell;
    cell.nat = 1;
    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.sccs_config.cavity.epsilon_bulk = 78.3;

    surchem solvent;
    solvent.set_parameters(parameters);
    ModuleBase::matrix unused_vloc;
    ModuleBase::matrix force(1, 3);
    EXPECT_THROW(solvent.cal_force_sol(cell, nullptr, unused_vloc, 1, force),
                 std::invalid_argument);
}

TEST(SolForce, ConvertsPointIonPcc2dForceFromHartreeToRydberg)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 2.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    const double volume = 2.0 * length * length * length;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    basis.collect_uniqgg();

    UnitCell cell;
    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = volume;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 2;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 2;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 0.5;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.4, 0.4, 0.5));
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.6, 0.6, 0.5));

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 1.0;
    parameters.expected_ionic_charge = 1.0;
    parameters.normalization_tolerance = 1.0e-10;
    parameters.sccs_config.cavity.density_min = 1.0e-2;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 1.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.sccs_config.max_iterations = 20;
    parameters.sccs_config.mixing = 0.7;
    parameters.sccs_config.tolerance_rms = 1.0e-14;
    parameters.sccs_config.tolerance_max = 1.0e-14;

    surchem solvent;
    solvent.set_parameters(parameters);
    std::vector<double> electron_density(basis.nrxx, 1.0 / volume);
    const double* density_channels[1] = {electron_density.data()};
    std::vector<double> local_potential(basis.nrxx, 0.0);
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);

    ModuleBase::matrix unused_vloc(1, basis.ngg);
    ModuleBase::matrix force(2, 3);
    solvent.cal_force_sol(cell, &basis, unused_vloc, 1, force);

    ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(lattice, length, 1.0e-10);
    const std::vector<double> positions_y{4.0, 6.0};
    const std::vector<double> masses{1.0, 1.0};
    geometry.origin_y = ModuleSccs::pcc_2d_system_center_y(
        positions_y,
        masses,
        geometry.parameters.cell_length_y);
    for (int atom = 0; atom < 2; ++atom)
    {
        ModuleSccs::PointCharge point;
        point.charge = 0.5;
        point.position = cell.atoms[0].tau[atom] * length;
        const ModuleBase::Vector3<double> expected_hartree
            = ModuleSccs::pcc_2d_point_charge_force(
                solvent.sccs_result().point_solute_moments_2d,
                point,
                geometry);
        EXPECT_NEAR(force(atom, 0), 2.0 * expected_hartree.x, 1.0e-14);
        EXPECT_NEAR(force(atom, 1), 2.0 * expected_hartree.y, 1.0e-14);
        EXPECT_NEAR(force(atom, 2), 2.0 * expected_hartree.z, 1.0e-14);
    }
}

TEST(SolForce, MatchesFixedElectronDensityReactionEnergyDerivative)
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
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.38, 0.47, 0.51));

    ModuleBase::matrix radial_local_potential(1, basis.ngg);
    for (int radial_index = 0; radial_index < basis.ngg; ++radial_index)
    {
        const double gg = basis.gg_uniq[radial_index];
        if (gg == 0.0)
        {
            continue;
        }
        const double coulomb_rydberg
            = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * gg);
        radial_local_potential(0, radial_index)
            = -coulomb_rydberg / cell.omega * std::exp(-0.3 * cell.tpiba2 * gg);
    }

    std::vector<double> electron_density(basis.nrxx);
    double electron_count = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        const double dx = (static_cast<double>(ix) + 0.5) / basis.nx - 0.62;
        const double dy = (static_cast<double>(iy) + 0.5) / basis.ny - 0.48;
        const double dz = (static_cast<double>(iz) + 0.5) / basis.nz - 0.50;
        electron_density[ir]
            = std::exp(-35.0 * (dx * dx + dy * dy + dz * dz));
        electron_count += electron_density[ir];
    }
    electron_count *= cell.omega / static_cast<double>(basis.nxyz);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electron_density[ir] /= electron_count;
    }

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 1.0;
    parameters.expected_ionic_charge = 1.0;
    parameters.normalization_tolerance = 1.0e-10;
    parameters.sccs_config.cavity.density_min = 1.0e-3;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 5.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Periodic;
    parameters.pcc_boundary = ModuleSccs::Boundary::Periodic;
    parameters.sccs_config.max_iterations = 300;
    parameters.sccs_config.mixing = 0.5;
    parameters.sccs_config.tolerance_rms = 1.0e-12;
    parameters.sccs_config.tolerance_max = 1.0e-10;

    surchem solvent;
    solvent.set_parameters(parameters);
    evaluate_reaction_energy(solvent,
                             cell,
                             basis,
                             radial_local_potential,
                             electron_density);
    ModuleBase::matrix force(1, 3);
    solvent.cal_force_sol(cell, &basis, radial_local_potential, 1, force);

    const double displacement = 1.0e-3;
    cell.atoms[0].tau[0].x += displacement / length;
    const double energy_plus
        = evaluate_reaction_energy(solvent,
                                   cell,
                                   basis,
                                   radial_local_potential,
                                   electron_density);
    cell.atoms[0].tau[0].x -= 2.0 * displacement / length;
    const double energy_minus
        = evaluate_reaction_energy(solvent,
                                   cell,
                                   basis,
                                   radial_local_potential,
                                   electron_density);
    const double finite_difference_force
        = -(energy_plus - energy_minus) / (2.0 * displacement);

    EXPECT_NEAR(0.5 * force(0, 0), finite_difference_force, 1.0e-4);
}

TEST(SolForce, NeutralAndChargedPcc2dMatchFixedDensityTotalEnergyDerivativeInXyz)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 2.0, 0.0,
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
    cell.omega = 2.0 * length * length * length;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.38, 0.47, 0.51));

    ModuleBase::matrix radial_local_potential(1, basis.ngg);
    for (int radial_index = 0; radial_index < basis.ngg; ++radial_index)
    {
        const double gg = basis.gg_uniq[radial_index];
        if (gg == 0.0)
        {
            continue;
        }
        const double coulomb_rydberg
            = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * gg);
        radial_local_potential(0, radial_index)
            = -coulomb_rydberg / cell.omega * std::exp(-0.3 * cell.tpiba2 * gg);
    }

    std::vector<double> electron_density(basis.nrxx);
    double electron_count = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const int iy = ir / basis.nplane - ix * basis.ny;
        const int iz = ir % basis.nplane + basis.startz_current;
        const double dx = (static_cast<double>(ix) + 0.5) / basis.nx - 0.62;
        const double dy = (static_cast<double>(iy) + 0.5) / basis.ny - 0.43;
        const double dz = (static_cast<double>(iz) + 0.5) / basis.nz - 0.55;
        electron_density[ir]
            = std::exp(-35.0 * (dx * dx + dy * dy + dz * dz));
        electron_count += electron_density[ir];
    }
    electron_count *= cell.omega / static_cast<double>(basis.nxyz);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electron_density[ir] /= electron_count;
    }

    const double displacement = 1.0e-3;
    const double coordinate_scale[3] = {length, length, length};
    double* coordinates[3] = {&cell.atoms[0].tau[0].x,
                              &cell.atoms[0].tau[0].y,
                              &cell.atoms[0].tau[0].z};
    const double electron_counts[2] = {1.0, 0.8};
    for (int charge_case = 0; charge_case < 2; ++charge_case)
    {
        std::vector<double> case_density(electron_density);
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            case_density[ir] *= electron_counts[charge_case];
        }

        SurchemParameters parameters;
        parameters.use_sccs = true;
        parameters.expected_electron_count = electron_counts[charge_case];
        parameters.expected_ionic_charge = 1.0;
        parameters.normalization_tolerance = 1.0e-10;
        parameters.sccs_config.cavity.density_min = 1.0e-3;
        parameters.sccs_config.cavity.density_max = 2.0e-2;
        parameters.sccs_config.cavity.epsilon_bulk = 1.1;
        parameters.sccs_config.surface_regularization = 1.0e-6;
        parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc2d;
        parameters.sccs_config.max_iterations = 500;
        parameters.sccs_config.mixing = 0.5;
        parameters.sccs_config.tolerance_rms = 1.0e-12;
        parameters.sccs_config.tolerance_max = 1.0e-10;

        surchem solvent;
        solvent.set_parameters(parameters);
        evaluate_total_electrostatic_energy(solvent,
                                            cell,
                                            basis,
                                            radial_local_potential,
                                            case_density);
        ModuleBase::matrix force(1, 3);
        solvent.cal_force_sol(cell, &basis, radial_local_potential, 1, force);

        for (int direction = 0; direction < 3; ++direction)
        {
            *coordinates[direction] += displacement / coordinate_scale[direction];
            const double energy_plus
                = evaluate_total_electrostatic_energy(solvent,
                                                      cell,
                                                      basis,
                                                      radial_local_potential,
                                                      case_density);
            *coordinates[direction] -= 2.0 * displacement / coordinate_scale[direction];
            const double energy_minus
                = evaluate_total_electrostatic_energy(solvent,
                                                      cell,
                                                      basis,
                                                      radial_local_potential,
                                                      case_density);
            *coordinates[direction] += displacement / coordinate_scale[direction];
            const double finite_difference_force
                = -(energy_plus - energy_minus) / (2.0 * displacement);
            EXPECT_NEAR(0.5 * force(0, direction), finite_difference_force, 2.0e-4)
                << "charge case " << charge_case << " direction " << direction;
        }
    }
}

TEST(SolForce, RejectsIncorrectOutputShape)
{
    UnitCell cell;
    cell.nat = 2;
    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.sccs_config.cavity.epsilon_bulk = 1.0;

    surchem solvent;
    solvent.set_parameters(parameters);
    ModuleBase::matrix unused_vloc;
    ModuleBase::matrix force(1, 3);
    EXPECT_THROW(solvent.cal_force_sol(cell, nullptr, unused_vloc, 1, force),
                 std::invalid_argument);
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
