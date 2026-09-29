#include "../pcc/pcc_boundary.h"

#include "gtest/gtest.h"

#include <stdexcept>

namespace
{

TEST(PccBoundary, ParsesNamesCaseInsensitively)
{
    EXPECT_EQ(ModulePcc::parse_boundary("periodic"), ModulePcc::Boundary::Periodic);
    EXPECT_EQ(ModulePcc::parse_boundary("PCC_0D"), ModulePcc::Boundary::Pcc0d);
    EXPECT_EQ(ModulePcc::parse_boundary("PCC-2D"), ModulePcc::Boundary::Pcc2d);
    EXPECT_THROW(ModulePcc::parse_boundary("pcc_1d"), std::invalid_argument);
}

} // namespace
