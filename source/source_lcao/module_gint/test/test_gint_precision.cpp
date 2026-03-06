#include "gtest/gtest.h"

#include "../gint_precision.h"

TEST(GintPrecisionTest, ScopedExecConfigRestoresPreviousPrecision)
{
    EXPECT_EQ(ModuleGint::current_exec_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);

    {
        const ModuleGint::ScopedExecConfig fp32_scope(
            ModuleGint::GintExecConfig{ModuleGint::GintRealPrecision::fp32});
        EXPECT_EQ(ModuleGint::current_exec_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);

        {
            const ModuleGint::ScopedExecConfig fp64_scope(
                ModuleGint::GintExecConfig{ModuleGint::GintRealPrecision::fp64});
            EXPECT_EQ(ModuleGint::current_exec_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);
        }

        EXPECT_EQ(ModuleGint::current_exec_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);
    }

    EXPECT_EQ(ModuleGint::current_exec_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);
}
