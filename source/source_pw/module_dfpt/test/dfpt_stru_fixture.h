#ifndef DFPT_STRU_FIXTURE_H
#define DFPT_STRU_FIXTURE_H

#include "source_cell/unitcell.h"

#include "gtest/gtest.h"
#include <string>
#include <vector>

// Shared gtest fixture for building a minimal cubic UnitCell from a
// hand-written structure table (abbreviated from
// module_symmetry/test/symm_test.cpp and klist_test.cpp). Used by the
// MPI-side DFPT tests that drive the QList / DFPT_PW wiring
// (dfpt_pw_data_test.cpp, dfpt_pw_run_test.cpp).
//
// All members touched by construct_ucell (UnitCell geometry fields and
// the Atom label/na/tau/taud vectors) are public, so tests include the
// cell headers normally.

struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group;    // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

class DFPTStruTestFixture : public testing::Test
{
  protected:
    UnitCell ucell;
    std::vector<stru_> stru_lib;

    DFPTStruTestFixture();

    void construct_ucell(stru_& stru);
    void ClearUcell();
};

#endif // DFPT_STRU_FIXTURE_H
