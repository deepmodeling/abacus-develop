#ifndef DFPT_STRU_FIXTURE_H
#define DFPT_STRU_FIXTURE_H

#include <string>
#include <vector>
#include "source_cell/unitcell.h"
#include "gtest/gtest.h"

// Shared gtest fixture for building a minimal cubic UnitCell from a
// hand-written structure table (abbreviated from
// module_symmetry/test/symm_test.cpp and klist_test.cpp). Used by the
// MPI-side DFPT tests that drive the QList / DFPT_PW wiring
// (dfpt_pw_data_test.cpp, dfpt_pw_run_test.cpp).
//
// NOTE ON INCLUDE ORDER: every test that needs UnitCell private members
// includes the cell headers with `#define private public` BEFORE this
// header; the include guards then keep the fixture header's own includes
// inert. The fixture implementation (dfpt_stru_fixture.cpp) only touches
// public members, so it compiles without the define.

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
