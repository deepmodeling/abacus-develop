#ifndef SCCS_PW_CHARGE_H
#define SCCS_PW_CHARGE_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <vector>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleSccs
{

double validate_cubic_cell(const ModuleBase::Matrix3& lattice_vectors,
                           double lattice_scale,
                           double relative_tolerance);

ModuleBase::Vector3<double> cell_center(const ModuleBase::Matrix3& lattice_vectors,
                                       double lattice_scale);

std::vector<ModuleBase::Vector3<double>> pw_grid_positions(
    const ModulePW::PW_Basis& basis,
    const ModuleBase::Matrix3& lattice_vectors,
    double lattice_scale);

std::vector<double> ionic_charge_from_local_potential(
    const std::vector<double>& local_potential_rydberg,
    double total_ionic_charge,
    double cell_volume,
    double tpiba,
    const ModulePW::PW_Basis& basis);

} // namespace ModuleSccs

#endif
