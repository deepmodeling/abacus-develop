/**
 * @file cell_tools.h
 * @brief Free function tools for extracting cell/atom information.
 */
#ifndef CELL_TOOLS_H
#define CELL_TOOLS_H

#include <string>
#include <vector>

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
}

#endif // CELL_TOOLS_H
