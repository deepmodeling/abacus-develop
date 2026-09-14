#include "gtest/gtest.h"

#include "../gint_info.h"

TEST(GintPrecisionTest, GintInfoStoresExecPrecision)
{
    UnitCell ucell;
    ModuleGint::GintInfo* info = ModuleGint::GintInfo::make_test_instance_ptr(ucell, {});

    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->get_exec_precision(), ModuleGint::GintPrecision::fp64);

    info->set_exec_precision(ModuleGint::GintPrecision::fp32);
    EXPECT_EQ(info->get_exec_precision(), ModuleGint::GintPrecision::fp32);
}
