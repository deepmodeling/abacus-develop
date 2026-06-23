#ifndef SYNTHETIC_NEIGHBOR_UNITCELL_H
#define SYNTHETIC_NEIGHBOR_UNITCELL_H

#include "prepare_unitcell.h"

#include <string>
#include <valarray>
#include <vector>

struct SyntheticNeighborCase
{
    std::string name;
    UcellTestPrepare prepare;
    std::vector<double> radii;
};

inline std::vector<SyntheticNeighborCase> make_synthetic_neighbor_cases()
{
    std::vector<SyntheticNeighborCase> cases;

    // Orthogonal single-element cell with several atoms distributed across
    // different boxes. It is small enough for cheap full pair comparison.
    cases.push_back(SyntheticNeighborCase{
        "cubic_multi_atom",
        UcellTestPrepare("sc",
                         2,
                         false,
                         false,
                         false,
                         "None",
                         10.0,
                         std::valarray<double>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                         std::vector<std::string>{"Si"},
                         std::vector<std::string>{"Si.upf"},
                         std::vector<std::string>{"auto"},
                         std::vector<std::string>{"Si.orb"},
                         std::valarray<int>{4},
                         std::vector<double>{28.085},
                         "Cartesian",
                         std::valarray<double>{0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.9, 0.0, 0.9, 0.9, 0.9}),
        std::vector<double>{0.95, 1.6}});

    // Skewed triclinic cell exercises direct/cartesian conversion and periodic
    // image shifts where box coordinates are not axis-aligned with the lattice.
    cases.push_back(SyntheticNeighborCase{
        "skewed_triclinic",
        UcellTestPrepare("triclinic",
                         2,
                         false,
                         false,
                         false,
                         "None",
                         10.0,
                         std::valarray<double>{1.0, 0.0, 0.0, 0.28, 0.96, 0.0, 0.16, 0.18, 1.08},
                         std::vector<std::string>{"Si"},
                         std::vector<std::string>{"Si.upf"},
                         std::vector<std::string>{"auto"},
                         std::vector<std::string>{"Si.orb"},
                         std::valarray<int>{4},
                         std::vector<double>{28.085},
                         "Cartesian",
                         std::valarray<double>{0.0, 0.0, 0.0, 0.7, 0.2, 0.0, 0.1, 0.8, 0.4, 0.9, 0.9, 0.9}),
        std::vector<double>{0.85, 1.5}});

    // Multi-element cell verifies that pair keys preserve type/natom identity
    // and not just geometric distance.
    cases.push_back(SyntheticNeighborCase{
        "multi_element",
        UcellTestPrepare("sc",
                         2,
                         false,
                         false,
                         false,
                         "None",
                         10.0,
                         std::valarray<double>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                         std::vector<std::string>{"Si", "O"},
                         std::vector<std::string>{"Si.upf", "O.upf"},
                         std::vector<std::string>{"auto", "auto"},
                         std::vector<std::string>{"Si.orb", "O.orb"},
                         std::valarray<int>{2, 3},
                         std::vector<double>{28.085, 15.999},
                         "Cartesian",
                         std::valarray<double>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
                                               0.0, 1.0, 1.0, 0.2, 0.5, 0.5, 1.0}),
        std::vector<double>{1.05, 1.8}});

    return cases;
}

#endif
