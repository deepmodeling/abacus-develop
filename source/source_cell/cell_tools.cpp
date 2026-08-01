/**
 * @file cell_tools.cpp
 * @brief Implementation of cell tool free functions.
 */
#include "cell_tools.h"

namespace unitcell
{
    std::vector<std::string> get_atomLabels(const Atom* atoms, const int ntype)
    {
        std::vector<std::string> atomLabels(ntype);
        for (int it = 0; it < ntype; it++)
        {
            atomLabels[it] = atoms[it].label;
        }
        return atomLabels;
    }

    std::vector<int> get_atomCounts(const Atom* atoms, const int ntype)
    {
        std::vector<int> atomCounts(ntype);
        for (int it = 0; it < ntype; it++)
        {
            atomCounts[it] = atoms[it].na;
        }
        return atomCounts;
    }

    std::vector<std::vector<int>> get_lnchiCounts(const Atom* atoms, const int ntype)
    {
        std::vector<std::vector<int>> lnchiCounts(ntype);
        for (int it = 0; it < ntype; it++)
        {
            lnchiCounts[it].resize(atoms[it].nwl + 1);
            for (int L = 0; L < atoms[it].nwl + 1; L++)
            {
                lnchiCounts[it][L] = atoms[it].l_nchi[L];
            }
        }
        return lnchiCounts;
    }
}
