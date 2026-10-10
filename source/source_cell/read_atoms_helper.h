#ifndef READ_ATOMS_HELPER_H
#define READ_ATOMS_HELPER_H

#include <fstream>
#include <string>
#include "unitcell.h"
#include "source_base/vector3.h"
#include "source_base/matrix3.h"

namespace unitcell {

/**
 * @brief Validate coordinate system type
 * @param Coordinate The coordinate system string to validate
 * @param ofs_warning Output stream for warnings
 * @return true if valid, false otherwise
 */
bool validate_coordinate_system(const std::string& Coordinate,
                                std::ofstream& ofs_warning);

/**
 * @brief Allocate and initialize atom property vectors
 * @param atom The atom object to allocate properties for
 * @param na Number of atoms
 * @param mass Atomic mass
 */
void allocate_atom_properties(Atom& atom, int na);

/**
 * @brief Set atom movement constraints based on fixed_atoms parameter
 * @param atom The atom object
 * @param ia Atom index
 * @param mv Movement vector (1=movable, 0=fixed)
 */
void set_atom_movement_flags(Atom& atom, int ia,
                            const ModuleBase::Vector3<int>& mv,
                            const bool fixed_atoms);

/**
 * @brief Set default magnetization if not explicitly specified
 * @param ucell Unit cell object
 * @param nspin Number of spin components
 * @param ofs_running Output stream for running information
 */
void autoset_magnetization(UnitCell& ucell, int nspin,
                          std::ofstream& ofs_running);

/**
 * @brief Check whether every atom carries a zero collinear magnetic moment
 * @param ucell Unit cell object
 * @return true if |mag| of every atom is below the magnetization threshold
 */
bool is_magnetization_all_zero(const UnitCell& ucell);

/**
 * @brief Warn when a nspin=4 run starts from an all-zero magnetization vector
 *
 * The check uses m_loc_, the vector actually consumed by cal_ux(), and not
 * mag: the collinear projection applied when noncolin=0 can leave mag
 * non-zero while m_loc_ is already zero, so keying the check on mag would
 * miss exactly the silent zero-moment start this warning reports.
 *
 * @param ucell Unit cell object
 * @param ofs_running Output stream for running information
 */
void warn_zero_magnetization_nspin4(const UnitCell& ucell,
                                    std::ofstream& ofs_running);

/**
 * @brief Report, once per atom type, that x/y magnetization was discarded
 * @param label Element label of the atom type
 * @param n_discarded Number of atoms of this type whose x/y moment was dropped
 * @param na Number of atoms of this type
 * @param ofs_running Output stream for running information
 */
void warn_xy_magnetization_ignored(const std::string& label,
                                   int n_discarded,
                                   int na,
                                   std::ofstream& ofs_running);

/**
 * @brief Perform final validation and output
 * @param ucell Unit cell object
 * @param ofs_running Output stream for running information
 * @param ofs_warning Output stream for warnings
 * @return true if validation passes, false otherwise
 */
bool finalize_atom_positions(UnitCell& ucell,
                            std::ofstream& ofs_running,
                            std::ofstream& ofs_warning,
                            const std::string& calculation,
                            const std::string& esolver_type);

/**
 * @brief Calculate lattice center for different centering modes
 * @param latvec Lattice vectors
 * @param center_mode Centering mode: "xy", "xz", "yz", or "xyz"
 * @return Lattice center coordinates
 */
ModuleBase::Vector3<double> calculate_lattice_center(
    const ModuleBase::Matrix3& latvec,
    const std::string& center_mode);

/**
 * @brief Convert between different coordinate systems
 * @param atom The atom object
 * @param ia Atom index
 * @param Coordinate Coordinate system type
 * @param v Input position vector
 * @param latvec Lattice vectors
 * @param lat0 Lattice constant
 * @param latcenter Lattice center (output parameter)
 */
void transform_atom_coordinates(Atom& atom, int ia,
                               const std::string& Coordinate,
                               const ModuleBase::Vector3<double>& v,
                               const ModuleBase::Matrix3& latvec,
                               double lat0,
                               ModuleBase::Vector3<double>& latcenter);

/**
 * @brief Convert between magnetization representations and output
 * @param atom The atom object
 * @param it Atom type index
 * @param ia Atom index
 * @param nspin Number of spin components
 * @param input_vec_mag Whether vector magnetization was input
 * @param input_angle_mag Whether angle magnetization was input
 * @param ofs_running Output stream for running information
 * @return true if the x/y components were discarded by the collinear
 *         projection applied when nspin=4 and noncolin=0
 */
bool process_magnetization(Atom& atom, int it, int ia,
                          int nspin, bool input_vec_mag,
                          bool input_angle_mag,
                          std::ofstream& ofs_running,
                          const bool noncolin);

/**
 * @brief Parse optional atom properties (mag, angle1, angle2, lambda, sc, m, v)
 * @param ifpos Input file stream
 * @param atom The atom object
 * @param ia Atom index
 * @param mv Movement vector (output parameter)
 * @param input_vec_mag Whether vector magnetization was input (output parameter)
 * @param input_angle_mag Whether angle magnetization was input (output parameter)
 * @param set_element_mag_zero Whether to reset element magnetization (output parameter)
 * @return true if parsing succeeds, false otherwise
 */
bool parse_atom_properties(std::ifstream& ifpos,
                          Atom& atom, int ia,
                          ModuleBase::Vector3<int>& mv,
                          bool& input_vec_mag,
                          bool& input_angle_mag,
                          bool& set_element_mag_zero);

/**
 * @brief Read atom type metadata (label, magnetization, orbital info, atom count)
 * @param it Atom type index
 * @param ucell Unit cell object
 * @param ifpos Input file stream
 * @param ofs_running Output stream for running information
 * @param ofs_warning Output stream for warnings
 * @param set_element_mag_zero Whether to reset element magnetization (output parameter)
 * @return true if reading succeeds, false otherwise
 */
bool read_atom_type_header(int it, UnitCell& ucell,
                          std::ifstream& ifpos,
                          std::ofstream& ofs_running,
                          std::ofstream& ofs_warning,
                          bool& set_element_mag_zero,
                          const std::string& basis_type,
                          const std::string& orbital_dir,
                          const std::string& init_wfc,
                          const double onsite_radius);

} // namespace unitcell

#endif // READ_ATOMS_HELPER_H
