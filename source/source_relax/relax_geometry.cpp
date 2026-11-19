#include "relax_geometry.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_common.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_cell/unitcell.h"
#include "source_cell/update_cell.h"

#include <cassert>
#include <cmath>

namespace relax_geometry
{
double project_ionic_forces(const UnitCell &ucell, const ModuleBase::matrix &force, ModuleBase::matrix &projected_force)
{
    const double to_ev_angstrom = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;
    projected_force = force * to_ev_angstrom;
    int iat = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        const Atom &atom = ucell.atoms[it];
        for (int ia = 0; ia < atom.na; ++ia, ++iat)
        {
            for (int axis = 0; axis < 3; ++axis)
            {
                if (atom.mbl[ia][axis] != 1)
                {
                    projected_force(iat, axis) = 0.0;
                }
            }
        }
    }
    assert(iat == force.nr);
    return projected_force.absmax();
}

double project_cell_stress(
    const UnitCell &ucell,
    const ModuleBase::matrix &stress,
    const std::string &fixed_axes,
    ModuleBase::matrix &projected_virial)
{
    const double virial_scale = ucell.omega * ModuleBase::Ry_to_eV;
    projected_virial = stress * virial_scale;
    if (fixed_axes == "shape")
    {
        const double mean = projected_virial.trace_on() / 3.0;
        projected_virial.zero_out();
        for (int axis = 0; axis < 3; ++axis)
        {
            projected_virial(axis, axis) = mean;
        }
    }
    else if (fixed_axes == "volume")
    {
        const double mean = projected_virial.trace_on() / 3.0;
        for (int axis = 0; axis < 3; ++axis)
        {
            projected_virial(axis, axis) -= mean;
        }
    }
    else if (fixed_axes != "None")
    {
        // A locked axis fixes a lattice vector, i.e. a row of latvec.
        ModuleBase::matrix lattice_virial = ucell.latvec.to_matrix() * projected_virial;
        for (int row = 0; row < 3; ++row)
        {
            if (ucell.lat_axis_free[row] == 0)
            {
                for (int col = 0; col < 3; ++col)
                {
                    lattice_virial(row, col) = 0.0;
                }
            }
        }
        projected_virial = ucell.GT.to_matrix() * lattice_virial;
    }

    const double to_kbar = ModuleBase::RYDBERG_SI / std::pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    return projected_virial.absmax() / virial_scale * to_kbar;
}

void apply_lattice_step(
    UnitCell &ucell,
    ModuleBase::Matrix3 displacement,
    bool fixed_volume,
    bool fixed_ibrav)
{
    // Remove any residual displacement of a locked lattice vector.
    if (ucell.lat_axis_free[0] == 0)
    {
        displacement.e11 = displacement.e12 = displacement.e13 = 0.0;
    }
    if (ucell.lat_axis_free[1] == 0)
    {
        displacement.e21 = displacement.e22 = displacement.e23 = 0.0;
    }
    if (ucell.lat_axis_free[2] == 0)
    {
        displacement.e31 = displacement.e32 = displacement.e33 = 0.0;
    }
    ucell.latvec += displacement;

    if (fixed_volume)
    {
        const double volume = std::abs(ucell.latvec.Det()) * std::pow(ucell.lat0, 3);
        ucell.latvec *= std::pow(ucell.omega / volume, 1.0 / 3.0);
    }
    if (fixed_ibrav)
    {
        unitcell::remake_cell(ucell.lat);
    }
}

void apply_ionic_step(
    UnitCell &ucell,
    const ModuleBase::Matrix3 &old_gt,
    ModuleBase::matrix displacement,
    double fraction)
{
    assert(displacement.nr == ucell.nat && displacement.nc == 3);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        ModuleBase::Vector3<double> cart;
        for (int axis = 0; axis < 3; ++axis)
        {
            cart[axis] = displacement(iat, axis) / ModuleBase::BOHR_TO_A / ucell.lat0;
        }
        const ModuleBase::Vector3<double> direct = cart * old_gt;
        const auto &movable = ucell.atoms[ucell.iat2it[iat]].mbl[ucell.iat2ia[iat]];
        for (int axis = 0; axis < 3; ++axis)
        {
            displacement(iat, axis) = movable[axis] == 1 ? direct[axis] * fraction : 0.0;
        }
    }

    if (ModuleSymmetry::Symmetry::symm_flag && ucell.symm.all_mbl && ucell.symm.nrotk > 0)
    {
        ucell.symm.symmetrize_vec3_nat(displacement.c);
    }
    unitcell::update_pos_taud(ucell.lat, displacement.c, ucell.ntype, ucell.nat, ucell.atoms);
}

void refresh_cell(UnitCell &ucell, int nspin, std::ofstream &ofs_running)
{
    ucell.a1.x = ucell.latvec.e11;
    ucell.a1.y = ucell.latvec.e12;
    ucell.a1.z = ucell.latvec.e13;
    ucell.a2.x = ucell.latvec.e21;
    ucell.a2.y = ucell.latvec.e22;
    ucell.a2.z = ucell.latvec.e23;
    ucell.a3.x = ucell.latvec.e31;
    ucell.a3.y = ucell.latvec.e32;
    ucell.a3.z = ucell.latvec.e33;

#ifdef __MPI
    // distribute lattice vectors.
    Parallel_Common::bcast_double(ucell.latvec.e11);
    Parallel_Common::bcast_double(ucell.latvec.e12);
    Parallel_Common::bcast_double(ucell.latvec.e13);
    Parallel_Common::bcast_double(ucell.latvec.e21);
    Parallel_Common::bcast_double(ucell.latvec.e22);
    Parallel_Common::bcast_double(ucell.latvec.e23);
    Parallel_Common::bcast_double(ucell.latvec.e31);
    Parallel_Common::bcast_double(ucell.latvec.e32);
    Parallel_Common::bcast_double(ucell.latvec.e33);

    // distribute lattice vectors.
    Parallel_Common::bcast_double(ucell.a1.x);
    Parallel_Common::bcast_double(ucell.a1.y);
    Parallel_Common::bcast_double(ucell.a1.z);
    Parallel_Common::bcast_double(ucell.a2.x);
    Parallel_Common::bcast_double(ucell.a2.y);
    Parallel_Common::bcast_double(ucell.a2.z);
    Parallel_Common::bcast_double(ucell.a3.x);
    Parallel_Common::bcast_double(ucell.a3.y);
    Parallel_Common::bcast_double(ucell.a3.z);
#endif

    ucell.omega = std::abs(ucell.latvec.Det()) * ucell.lat0 * ucell.lat0 * ucell.lat0;

    ucell.GT = ucell.latvec.Inverse();
    ucell.G = ucell.GT.Transpose();
    ucell.GGT = ucell.G * ucell.GT;
    ucell.invGGT = ucell.GGT.Inverse();

    // prepare something for next SCF
    unitcell::setup_cell_after_vc(ucell, ofs_running, nspin);
    ModuleBase::GlobalFunc::DONE(ofs_running, "SETUP UNITCELL");
}
} // namespace relax_geometry
