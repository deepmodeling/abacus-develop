/**
 * @file unitcell_stats.cpp
 * @brief Per-species counting helpers of UnitCell:
 *        get_atom_Counts, get_orbital_Counts, get_lnchi_Counts.
 */
#include "unitcell.h"

std::map<int, int> UnitCell::get_atom_Counts() const {
    std::map<int, int> atomCounts;
    for (int it = 0; it < this->ntype; it++) {
        atomCounts.insert(std::pair<int, int>(it, this->atoms[it].na));
    }
    return atomCounts;
}

std::map<int, int> UnitCell::get_orbital_Counts() const {
    std::map<int, int> orbitalCounts;
    for (int it = 0; it < this->ntype; it++) {
        orbitalCounts.insert(std::pair<int, int>(it, this->atoms[it].nw));
    }
    return orbitalCounts;
}

std::map<int, std::map<int, int>> UnitCell::get_lnchi_Counts() const {
    std::map<int, std::map<int, int>> lnchiCounts;
    for (int it = 0; it < this->ntype; it++) {
        for (int L = 0; L < this->atoms[it].nwl + 1; L++) {
            // Check if the key 'it' exists in the outer map
            if (lnchiCounts.find(it) == lnchiCounts.end()) {
                // If it doesn't exist, initialize an empty inner map
                lnchiCounts[it] = std::map<int, int>();
            }
            int l_nchi = this->atoms[it].l_nchi[L];
            // Insert the key-value pair into the inner map
            lnchiCounts[it].insert(std::pair<int, int>(L, l_nchi));
        }
    }
    return lnchiCounts;
}
