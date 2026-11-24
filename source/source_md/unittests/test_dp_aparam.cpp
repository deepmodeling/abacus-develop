#include "source_md/module_dp/dp_aparam.h"

#include "source_cell/module_neighlist/local_atom.h"

#include "gtest/gtest.h"

#include <stdexcept>

TEST(DPAparamTest, ReordersOwnedAndGhostByGlobalIdentity)
{
    std::vector<LocalAtom> owned(2);
    owned[0].type = 1;
    owned[0].type_index = 0;
    owned[1].type = 0;
    owned[1].type_index = 1;
    std::vector<LocalAtom> ghost(1);
    ghost[0].type = 0;
    ghost[0].type_index = 0;
    const std::vector<std::int64_t> counts{2, 1};
    const std::vector<double> input{10.0, 20.0, 30.0, 40.0, 50.0, 60.0};

    const std::vector<double> nloc = ModuleMD::DPDetail::arrange_aparam(input, owned, ghost, counts, 2, false);
    const std::vector<double> nall = ModuleMD::DPDetail::arrange_aparam(input, owned, ghost, counts, 2, true);
    EXPECT_EQ(nloc, (std::vector<double>{50.0, 60.0, 30.0, 40.0}));
    EXPECT_EQ(nall, (std::vector<double>{50.0, 60.0, 30.0, 40.0, 10.0, 20.0}));

    const std::vector<double> uniform{7.0, 8.0};
    const std::vector<double> repeated = ModuleMD::DPDetail::arrange_aparam(uniform, owned, ghost, counts, 2, true);
    EXPECT_EQ(repeated, (std::vector<double>{7.0, 8.0, 7.0, 8.0, 7.0, 8.0}));
}

TEST(DPAparamTest, ValidatesShapeAndIdentity)
{
    std::vector<LocalAtom> owned(1);
    owned[0].type = 0;
    owned[0].type_index = 0;
    const std::vector<LocalAtom> ghost;
    const std::vector<std::int64_t> counts{1};
    EXPECT_TRUE(ModuleMD::DPDetail::arrange_aparam({}, owned, ghost, counts, 0, false).empty());
    EXPECT_THROW(ModuleMD::DPDetail::arrange_aparam({1.0}, owned, ghost, counts, 0, false), std::invalid_argument);
    EXPECT_THROW(ModuleMD::DPDetail::arrange_aparam({1.0}, owned, ghost, counts, 2, false), std::invalid_argument);
    owned[0].type_index = 1;
    EXPECT_THROW(ModuleMD::DPDetail::arrange_aparam({1.0, 2.0}, owned, ghost, counts, 2, false), std::invalid_argument);
}
