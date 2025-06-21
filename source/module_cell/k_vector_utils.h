//
// Created by rhx on 25-6-3.
//

#ifndef K_VECTOR_UTILS_H
#define K_VECTOR_UTILS_H

#include "source_base/global_variable.h"
#include "source_base/matrix3.h"
#include "module_cell/unitcell.h"

class K_Vectors;

namespace KVectorUtils
{
void kvec_d2c(K_Vectors& kv, const ModuleBase::Matrix3& reciprocal_vec);

void kvec_c2d(K_Vectors& kv, const ModuleBase::Matrix3& latvec);

/**
 * @brief Sets both the direct and Cartesian k-vectors.
 *
 * This function sets both the direct and Cartesian k-vectors based on the input parameters.
 * It also checks the k-point type and sets the corresponding flags.
 *
 * @param kv The K_Vectors object containing the k-point information.
 * @param G The reciprocal lattice matrix.
 * @param R The real space lattice matrix.
 * @param skpt A string to store the k-point table.
 *
 * @return void
 *
 * @note If the k-point type is neither "Cartesian" nor "Direct", an error message will be printed.
 * @note The function sets the flags kd_done and kc_done to indicate whether the direct and Cartesian k-vectors have
 * been set, respectively.
 * @note The function also prints a table of the direct k-vectors and their weights.
 * @note If the function is called by the master process (MY_RANK == 0), the k-point table is also stored in the
 * string skpt.
 */
void set_both_kvec(K_Vectors& kv, const ModuleBase::Matrix3& G, const ModuleBase::Matrix3& R, std::string& skpt);

/**
 * @brief Sets up the k-points after a volume change.
 *
 * This function sets up the k-points after a volume change in the system.
 * It sets the Cartesian and direct k-vectors based on the new reciprocal and real space lattice vectors.
 *
 * @param kv The K_Vectors object containing the k-point information.
 * @param nspin_in The number of spins. 1 for non-spin-polarized calculations and 2 for spin-polarized calculations.
 * @param reciprocal_vec The new reciprocal lattice matrix.
 *
 * @return void
 *
 * @note The function first sets the number of spins (nspin) to the input value.
 * @note The direct k-vectors have been set (kd_done = true) but the Cartesian k-vectors have not (kc_done =
 * false) after a volume change. The function calculates the Cartesian k-vectors by multiplying the direct k-vectors
 * with the reciprocal lattice matrix.
 * @note The function also prints a table of the direct k-vectors and their weights.
 * @note The function calls the print_klists function to print the k-points in both Cartesian and direct
 * coordinates.
 */
void set_after_vc(K_Vectors& kv, const int& nspin, const ModuleBase::Matrix3& G);

/**
 * @brief Prints the k-points in both Cartesian and direct coordinates.
 *
 * This function prints the k-points in both Cartesian and direct coordinates to the output file stream.
 * The output includes the index, x, y, and z coordinates, and the weight of each k-point.
 *
 * @param ofs The output file stream to which the k-points are printed.
 *
 * @return void
 *
 * @note The function first checks if the total number of k-points (nkstot) is less than the number of k-points for
 * the current spin (nks). If so, it prints an error message and quits.
 * @note The function prints the k-points in a table format, with separate tables for Cartesian and direct
 * coordinates.
 * @note The function uses the FmtCore::format function to format the output.
 */
void print_klists(const K_Vectors& kv, std::ofstream& ofs);

// step 3 : mpi kpoints information.

/**
 * @brief Distributes k-points among MPI processes.
 *
 * This function distributes the k-points among the MPI processes. Each process gets a subset of the k-points to
 * work on. The function also broadcasts various variables related to the k-points to all processes.
 *
 * @param kv The K_Vectors object containing the k-point information.
 *
 * @return void
 *
 * @note This function is only compiled and used if MPI is enabled.
 * @note The function assumes that the number of k-points (nkstot) is greater than 0.
 * @note The function broadcasts the flags kc_done and kd_done, the number of spins (nspin), the total number of
 * k-points (nkstot), the full number of k-points (nkstot_full), the Monkhorst-Pack grid (nmp), the k-point offsets
 * (koffset), and the segment IDs of the k-points (kl_segids).
 * @note The function also broadcasts the indices of the k-points (isk), their weights (wk), and their Cartesian and
 * direct coordinates (kvec_c and kvec_d).
 * @note If a process has no k-points to work on, the function will quit with an error message.
 */
#ifdef __MPI
void kvec_mpi_k(K_Vectors& kv);
#endif // __MPI

/**
 * @brief Generates irreducible k-points in the Brillouin zone considering symmetry operations.
 *
 * This function calculates the irreducible k-points (IBZ) from the given k-points, taking into
 * account the symmetry of the unit cell. It updates the symmetry-matched k-points and generates
 * the corresponding weight for each k-point.
 *
 * @param symm The symmetry information of the system.
 * @param use_symm A flag indicating whether to use symmetry operations.
 * @param skpt A string to store the formatted k-points information.
 * @param ucell The unit cell of the crystal.
 * @param match A boolean flag that indicates if the results matches the real condition.
 */
void kvec_ibz_kpoint(K_Vectors& kv,
                     const ModuleSymmetry::Symmetry& symm,
                     bool use_symm,
                     std::string& skpt,
                     const UnitCell& ucell,
                     bool& match);
} // namespace KVectorUtils

#endif // K_VECTOR_UTILS_H
