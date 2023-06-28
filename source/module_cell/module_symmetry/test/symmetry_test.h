#pragma once
#include "module_base/mathzone.h"
#include "../symmetry.h"
#include "gtest/gtest.h"

#define DOUBLETHRESHOLD 1e-8

struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group; // Schoenflies symbol
    std::string point_group_hm; // Hermann–Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
    std::string coordtype; // caltesian or direct
};

class SymmetryTest : public testing::Test
{
protected:
    UnitCell ucell;
    std::ofstream ofs_running;
    void construct_ucell(stru_& stru);
    void ClearUcell();
};