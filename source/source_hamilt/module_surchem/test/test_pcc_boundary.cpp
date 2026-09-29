#include "../pcc/pcc_boundary.h"

#include "gtest/gtest.h"

#include <stdexcept>

namespace
{

TEST(PccBoundary, ParsesNamesCaseInsensitively)
{
    EXPECT_EQ(ModuleSccs::parse_boundary("periodic"), ModuleSccs::Boundary::Periodic);
    EXPECT_EQ(ModuleSccs::parse_boundary("PCC_0D"), ModuleSccs::Boundary::Pcc0d);
    EXPECT_EQ(ModuleSccs::parse_boundary("PCC-2D"), ModuleSccs::Boundary::Pcc2d);
    EXPECT_THROW(ModuleSccs::parse_boundary("pcc_1d"), std::invalid_argument);
}

} // namespace
