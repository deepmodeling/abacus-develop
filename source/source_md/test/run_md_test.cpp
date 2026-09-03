#include "gtest/gtest.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private

#include "md_test_fixture.h"
#include "source_md/run_md.h"

TEST_F(MdTestBase, prepare_mdcell_from_unitcell)
{
    MDCell mdcell;
    Run_MD::prepare_mdcell(mdcell, ucell);

    EXPECT_EQ(mdcell.nat(), ucell.nat);
    EXPECT_EQ(mdcell.stru_file_metadata().species.size(), static_cast<std::size_t>(ucell.ntype));
}
