#include "source_lcao/module_deltaspin/deltaspin_state.h"

#include "gtest/gtest.h"

namespace
{

using ModuleBase::Vector3;
using spinconstrain::ScState;

TEST(ScStateCalEsconTest, CalEsconEmpty)
{
    ScState state;

    EXPECT_DOUBLE_EQ(state.cal_escon(), 0.0);
    EXPECT_DOUBLE_EQ(state.get_escon(), 0.0);
}

TEST(ScStateCalEsconTest, CalEsconSingleAtom)
{
    ScState state;
    state.set_atomCounts({{0, 1}});
    state.get_lambda().push_back(Vector3<double>(1.0, 1.0, 1.0));
    state.get_mi().push_back(Vector3<double>(1.0, 1.0, 1.0));

    EXPECT_DOUBLE_EQ(state.cal_escon(), -3.0);
    EXPECT_DOUBLE_EQ(state.get_escon(), -3.0);
}

TEST(ScStateCalEsconTest, CalEsconMultiAtom)
{
    ScState state;
    state.set_atomCounts({{0, 2}});
    state.get_lambda().push_back(Vector3<double>(1.0, 0.0, 0.0));
    state.get_lambda().push_back(Vector3<double>(0.0, 1.0, 0.0));
    state.get_mi().push_back(Vector3<double>(2.0, 0.0, 0.0));
    state.get_mi().push_back(Vector3<double>(0.0, 3.0, 0.0));

    EXPECT_DOUBLE_EQ(state.cal_escon(), -5.0);
    EXPECT_DOUBLE_EQ(state.get_escon(), -5.0);
}

} // namespace
