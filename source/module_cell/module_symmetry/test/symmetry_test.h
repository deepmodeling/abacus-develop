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
    std::vector<int> force_zero_iat;   // the index of atoms whose force should be zero
    std::map<int, int> force_oppo_iat;    // the index of atoms  pairs whose forces should be opposite
    std::vector<std::vector<int>> force_oppo_iat_xyz; //{ia1, ia2, xoppo(1)/eq(0), yoppo, zoppo}
    std::vector<std::pair<int, int>> stress_zero; //a set of elements in the stress tensor that should be zero
    std::vector<std::vector<std::pair<int, int>>> stress_eq; //a set of elements in the stress tensor that should be equal
};

class SymmetryTest : public testing::Test
{
protected:
    UnitCell ucell;
    std::ofstream ofs_running;
    void construct_ucell(stru_& stru);
    void ClearUcell();
};