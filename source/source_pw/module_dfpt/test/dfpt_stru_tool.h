#ifndef DFPT_STRU_TOOL_H
#define DFPT_STRU_TOOL_H

// Shared free-function utilities for the MPI-side DFPT tests
// (dfpt_pw_data_test / dfpt_pw_run_test): building a minimal UnitCell
// from a hand-written structure table (abbreviated from
// module_symmetry/test/symm_test.cpp and klist_test.cpp), in the same
// free-function test-tool pattern as test_serial/dfpt_test_tool.
// Every test file derives its own fixture from ::testing::Test, owns
// its UnitCell and calls these functions; there is no shared fixture
// base class.

#include <string>
#include <vector>

class UnitCell;

namespace dfpt_test
{

/// one atom species: label + cartesian coordinates of its sites
struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

/// hand-written structure table consumed by construct_ucell
struct stru_
{
    int ibrav;
    std::string point_group;    // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation
    std::string space_group;
    std::vector<double> cell;   // a1 | a2 | a3 (cartesian, row-major)
    std::vector<atomtype_> all_type;
};

/// The default test structure: a cubic 1x1x1 carbon cell (O_h / Pm-3m,
/// single atom at the origin).
stru_ cubic_carbon_stru();

/// Build the UnitCell from the structure table: lattice vectors and
/// latvec / G / GT, per-type tau / taud (cartesian -> direct through
/// Mathzone), and the iat2it / iat2ia running-index tables real cells
/// always carry. The caller owns the UnitCell (and its atoms array).
void construct_ucell(UnitCell& ucell, const stru_& stru);

} // namespace dfpt_test

#endif // DFPT_STRU_TOOL_H
