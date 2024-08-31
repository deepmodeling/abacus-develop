#include "module_cell/unitcell.h"
#include <string>
#include <vector>
#include <map>

namespace ModuleIO
{
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
 * For cif file from COD, which always contains experiments information and even SHELX
 * single-crystal refinement information, the CIF file is much more complicated.
 * see https://www.crystallography.net/cod/ for more information.
 * 
 * Design of this "class"
 * According to the CifParser implemented in pymatgen package.
 * 
 * Usage of this "class"
 * The usage of this class is simple, just call the static methods to_cif and from_cif.
 * There are also other utils such as vec_to_abc_angles, abc_angles_to_vec, etc.
 * A call similar with pymatgen implementation is also supported, in that case, a instance
 * of CifParser is bind with a CIF file, and the data can be accessed by the get() function.
 */
    class CifParser
    {
        public:
            CifParser() = delete; // I cannot understand why there should be a default constructor
            CifParser(const std::string& fcif) {}
            ~CifParser() {} // actually do not need to do anything...
            // some utils
            // general
            /**
             * @brief Convert the vector representation of lattice vectors to the cell A B C angles.
             * 
             * @param vec lattice vectors in the form of [a1, a2, a3, b1, b2, b3, c1, c2, c3]
             * @param abc_angles the output angles in the form of [a, b, c, alpha, beta, gamma]
             */
            static void vec_to_abc_angles(const double* vec, double* abc_angles);
            /**
             * @brief Convert the cell A B C angles to the vector representation of lattice vectors.
             * 
             * @param abc_angles the input angles in the form of [a, b, c, alpha, beta, gamma]
             * @param vec the output lattice vectors in the form of [a1, a2, a3, b1, b2, b3, c1, c2, c3]
             */
            static void abc_angles_to_vec(const double* abc_angles, double* vec);
            /**
             * @brief Convert the lattice vectors to the volume of the cell.
             * 
             * @param vec the input lattice vectors in the form of [a1, a2, a3, b1, b2, b3, c1, c2, c3]
             * @return double the volume of the cell
             */
            static double vec_to_volume(const double* vec);
            /**
             * @brief Convert the cell A B C angles to the volume of the cell.
             * 
             * @param abc_angles the input angles in the form of [a, b, c, alpha, beta, gamma]
             * @return double the volume of the cell
             */
            static double abc_angles_to_volume(const double* abc_angles);
            static double deg2rad(double deg) { return deg * M_PI / 180.0; }
            static double rad2deg(double rad) { return rad * 180.0 / M_PI; }
            /**
             * @brief Print the CIF file from the given information.
             * 
             * @param fcif the output cif file name
             * @param abc_angles cell A B C angles in the form of [a, b, c, alpha, beta, gamma]
             * @param natom the number of atoms
             * @param atom_site_labels the atom site labels in the form of [atom1, atom2, ...]
             * @param atom_site_fract_coords the fractional coordinates of the atoms in the form of [x1, y1, z1, x2, y2, z2, ...]
             * @param title the title of the CIF file
             * @param data_tag the data tag of the CIF file
             * @param atom_site_occups the occupancies of the atoms in the form of [occup1, occup2, ...]
             * @param cell_formula_units_z the cell formula units Z
             */
            static void to_cif(const std::string& fcif,
                               const double* abc_angles,
                               const int natom,
                               const std::string* atom_site_labels, // the one without numbers
                               const double* atom_site_fract_coords,
                               const std::string& title = "# generated by ABACUS",
                               const std::string& data_tag = "data_?",
                               const double* atom_site_occups = nullptr, // may be this will be useful after impementation of VCA?
                               const std::string& cell_formula_units_z = "1");
            // the version with both spacegroup symmetry and point group symmetry ready
            // not for now, because it is too complicated. However it is a walk-around
            // way to fix issue #4998
            // static void to_cif();
            /**
             * @brief Write CIF file with the whole UnitCell instance
             * 
             * @param fcif the output cif file name
             * @param ucell the input UnitCell instance
             * @param title the title of the CIF file
             * @param data_tag the data tag of the CIF file
             */
            static void to_cif(const std::string& fcif,
                               const UnitCell& ucell,
                               const std::string& title = "# generated by ABACUS",
                               const std::string& data_tag = "data_?");
            /**
             * @brief Read the CIF file and store the information in the map.
             * 
             * @param fcif the input cif file name
             * @param out the output map containing the information
             */
            static void from_cif(const std::string& fcif,
                                 std::map<std::string, std::vector<std::string>>& out);

            // not static :(
            std::vector<std::string> get(const std::string& key);
        // private:
            // not needed to be exposed
            static void _build_chem_formula(const int natom,
                                            const std::string* atom_site_labels,
                                            std::string& sum,
                                            std::string& structural);
            // interface to ABACUS UnitCell impl.
            static void _unpack_ucell(const UnitCell& ucell,    // because ucell is too heavy...
                                      std::vector<double>& veca,
                                      std::vector<double>& vecb,
                                      std::vector<double>& vecc,
                                      int& natom,
                                      std::vector<std::string>& atom_site_labels,
                                      std::vector<double>& atom_site_fract_coords);
            // split only those words out of specified enclose
            static std::vector<std::string> _split_outside_enclose(const std::string& in, 
                                                                   const std::string& delim,
                                                                   const std::vector<std::string>& enclose);
            static std::vector<std::string> _split_loop_block(const std::string& block);
            static std::map<std::string, std::vector<std::string>> _build_table(const std::vector<std::string>& keys,
                                                                                const std::vector<std::string>& values);
            static std::map<std::string, std::vector<std::string>> _build_block_data(const std::vector<std::string>& block);
        
        private:
            std::map<std::string, std::vector<std::string>> raw_;
    };
} // namespace ModuleIO
