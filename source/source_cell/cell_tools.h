/**
 * @file cell_tools.h
 * @brief Free function tools for extracting cell/atom information.
 */
#ifndef CELL_TOOLS_H
#define CELL_TOOLS_H

#include <string>
#include <vector>

#include "source_base/vector3.h"
#include "source_cell/atom_spec.h"

/**
 * @brief Free functions for extracting atom/orbital info from Atom array.
 */
namespace unitcell
{
    /// @brief Get atom labels for each atom type.
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @return vector of atom labels, one per type
    std::vector<std::string> get_atomLabels(const Atom* atoms, const int ntype);

    /// @brief Get atom counts (number of atoms) for each atom type.
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @return vector of atom counts, one per type
    std::vector<int> get_atomCounts(const Atom* atoms, const int ntype);

    /// @brief Get lnchi counts (number of chi functions per L) for each atom type.
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @return vector of lnchi counts, one vector per type
    std::vector<std::vector<int>> get_lnchiCounts(const Atom* atoms, const int ntype);

    /// @brief Get target magnetic moment for each atom (used by deltaspin).
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @param nat total number of atoms [in]
    /// @return vector of target magnetic moments, one per atom
    std::vector<ModuleBase::Vector3<double>> get_target_mag(const Atom* atoms,
                                                            const int ntype,
                                                            const int nat);

    /// @brief Get Lagrange multiplier for each atom (used by deltaspin).
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @param nat total number of atoms [in]
    /// @return vector of Lagrange multipliers, one per atom
    std::vector<ModuleBase::Vector3<double>> get_lambda(const Atom* atoms,
                                                         const int ntype,
                                                         const int nat);

    /// @brief Get constrain flag for each atom (used by deltaspin).
    /// @param atoms atom pointer [in]
    /// @param ntype number of atom types [in]
    /// @param nat total number of atoms [in]
    /// @return vector of constrain flags, one per atom
    std::vector<ModuleBase::Vector3<int>> get_constrain(const Atom* atoms,
                                                        const int ntype,
                                                        const int nat);
}

#endif // CELL_TOOLS_H
