#include "module_cell/unitcell.h"
#include <string>
#include <vector>

/**
 * @file cif_io.h
 * @brief This file contains the implementation of the functions for reading and writing CIF files.
 * 
 * Due to the space group and point group symmetry information are not always present, 
 * here will assume only P1 and C1 symmetry are present in the CIF file.
 * 
 * A typical CIF file with no symmetry (P1 and C1) is presented below, downloaded from
 * materials project: https://materialsproject.org/
 * 
 * # generated using pymatgen
 * data_C
 * _symmetry_space_group_name_H-M   'P 1'
 * _cell_length_a   2.46772428
 * _cell_length_b   2.46772428
 * _cell_length_c   8.68503800
 * _cell_angle_alpha   90.00000000
 * _cell_angle_beta   90.00000000
 * _cell_angle_gamma   120.00000758
 * _symmetry_Int_Tables_number   1
 * _chemical_formula_structural   C
 * _chemical_formula_sum   C4
 * _cell_volume   45.80317575
 * _cell_formula_units_Z   4
 * loop_
 *  _symmetry_equiv_pos_site_id
 *  _symmetry_equiv_pos_as_xyz
 *   1  'x, y, z'
 * loop_
 *  _atom_site_type_symbol
 *  _atom_site_label
 *  _atom_site_symmetry_multiplicity
 *  _atom_site_fract_x
 *  _atom_site_fract_y
 *  _atom_site_fract_z
 *  _atom_site_occupancy
 *   C  C0  1  0.00000000  0.00000000  0.75000000  1
 *   C  C1  1  0.00000000  0.00000000  0.25000000  1
 *   C  C2  1  0.33333300  0.66666700  0.75000000  1
 *   C  C3  1  0.66666700  0.33333300  0.25000000  1
 *
 */
namespace ModuleIO
{
    
} // namespace ModuleIO
