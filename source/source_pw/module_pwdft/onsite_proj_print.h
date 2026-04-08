#ifndef MODULEHAMILTPW_ONSITE_PROJ_PRINT_H
#define MODULEHAMILTPW_ONSITE_PROJ_PRINT_H

#include "source_cell/unitcell.h"
#include "source_base/global_variable.h"
#include <vector>
#include <string>
#include <complex>

namespace projectors {
namespace print {

/**
 * @brief Print orbital charge analysis
 * 
 * @param ucell Unit cell pointer
 * @param occs Occupation numbers
 * @param iat_nh Number of projectors per atom
 * @param atom_labels Atom labels
 */
void print_orb_chg(
    const UnitCell* ucell,
    const std::vector<std::complex<double>>& occs,
    const std::vector<int>& iat_nh,
    const std::vector<std::string>& atom_labels);

/**
 * @brief Print magnetism table
 * 
 * @param atom_labels Atom labels
 * @param mag_x Magnetic moment in x direction
 * @param mag_y Magnetic moment in y direction
 * @param mag_z Magnetic moment in z direction
 */
void print_mag_table(
    const std::vector<std::string>& atom_labels,
    const std::vector<double>& mag_x,
    const std::vector<double>& mag_y,
    const std::vector<double>& mag_z);

/**
 * @brief Print projector initialization status
 * 
 * @param it Atom type index
 * @param nproj_it Number of projectors for this type
 */
void print_proj_status(int it, int nproj_it);

} // namespace print
} // namespace projectors

#endif // MODULEHAMILTPW_ONSITE_PROJ_PRINT_H
