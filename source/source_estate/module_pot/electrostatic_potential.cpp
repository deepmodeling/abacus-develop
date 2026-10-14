#include "electrostatic_potential.h"
#include "h_hartree_pw.h"
#include "efield.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/module_parallel/para_world.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_io/module_parameter/input_parameter.h"

#include <cmath>
#include <stdexcept>

namespace elecstate
{
std::vector<double> electrostatic_potential(const UnitCell& cell, ModulePW::PW_Basis& basis,
                                           const Charge& charge, const double* fixed,
                                           const surchem& solvent, const Input_para& input,
                                           const std::vector<double>* pcc_potential)
{
    if (basis.nrxx > 0 && (!fixed || (input.imp_sol && !solvent.delta_phi)))
    { throw std::runtime_error("electrostatic reference requires initialized fixed and solvent potentials"); }
    const auto hartree = H_Hartree_pw::v_hartree(cell, &basis, input.nspin, charge.rho);
    ModuleBase::matrix dipole;
    const bool add_dipole = input.efield_flag && input.dip_cor_flag;
    if (add_dipole)
    { dipole = Efield::add_efield(cell, &basis, input.nspin, charge.rho, solvent); }
    std::vector<double> potential(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        potential[ir] = hartree(0, ir) + fixed[ir];
        if (add_dipole) { potential[ir] += dipole(0, ir); }
        if (pcc_potential) { potential[ir] += (*pcc_potential)[ir]; }
        if (input.imp_sol) { potential[ir] += solvent.delta_phi[ir]; }
    }
    return potential;
}

double planar_potential_average(const UnitCell& cell, const ModulePW::PW_Basis& basis,
                                const std::vector<double>& potential, int axis,
                                double minimum_A, double maximum_A)
{
    if (axis < 0 || axis > 2 || !std::isfinite(minimum_A) || !std::isfinite(maximum_A)
        || minimum_A < 0 || maximum_A <= minimum_A)
    { throw std::runtime_error("invalid electrostatic averaging axis/window"); }
    const double lattice[3][3] = {
        {cell.latvec.e11, cell.latvec.e12, cell.latvec.e13},
        {cell.latvec.e21, cell.latvec.e22, cell.latvec.e23},
        {cell.latvec.e31, cell.latvec.e32, cell.latvec.e33}};
    const int first = (axis + 1) % 3;
    const int second = (axis + 2) % 3;
    const double cross[3] = {
        lattice[first][1] * lattice[second][2] - lattice[first][2] * lattice[second][1],
        lattice[first][2] * lattice[second][0] - lattice[first][0] * lattice[second][2],
        lattice[first][0] * lattice[second][1] - lattice[first][1] * lattice[second][0]};
    double area2 = 0.0;
    double volume = 0.0;
    for (int i = 0; i < 3; ++i) { area2 += cross[i]*cross[i]; volume += lattice[axis][i]*cross[i]; }
    const double height_A = cell.lat0_angstrom * std::abs(volume) / std::sqrt(area2);
    if (!std::isfinite(height_A) || height_A <= 0.0 || maximum_A > height_A + 1.0e-10)
    { throw std::runtime_error("electrostatic averaging window lies outside the cell height"); }
    const int grid[3] = {basis.nx, basis.ny, basis.nz};
    if (grid[0] <= 0 || grid[1] <= 0 || grid[2] <= 0 || basis.nrxx < 0
        || (basis.nrxx > 0 && basis.nplane <= 0) || potential.size() != static_cast<std::size_t>(basis.nrxx))
    { throw std::runtime_error("invalid distributed FFT grid for electrostatic averaging"); }
    double reduction[3] = {0.0, 0.0, 0.0};
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int xy = ir / basis.nplane;
        const int indices[3] = {xy / basis.ny, xy % basis.ny, ir % basis.nplane + basis.startz_current};
        // FFT points and ABACUS cube origins are at i/n, not voxel centers.
        const double coordinate_A = double(indices[axis]) / grid[axis] * height_A;
        if (!std::isfinite(potential[ir])) { reduction[2] += 1.0; }
        if (coordinate_A >= minimum_A && coordinate_A < maximum_A)
        { reduction[0] += potential[ir]; reduction[1] += 1.0; }
    }
    // Reduce over this FFT grid's communicator without summing replicated pools.
#ifdef __MPI
    const auto grid_world = Parallel::ParaWorld::make_mpi("electrostatic_grid", basis.pool_world);
#else
    const auto grid_world = Parallel::ParaWorld::serial("electrostatic_grid");
#endif
    Parallel::reduce_all(reduction, 3, grid_world);
    if (reduction[2] || reduction[1] <= 0.0 || !std::isfinite(reduction[0]))
    { throw std::runtime_error("non-finite electrostatic potential or no FFT points in averaging window"); }
    return reduction[0] / reduction[1];
}
}
