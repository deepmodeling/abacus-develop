#include "gtest/gtest.h"

#include "../module_charge/gint_precision_controller.h"

TEST(GintPrecisionControllerTest, AutoModeSwitchesToFp64ImmediatelyWhenDrhoIsSmallEnough)
{
    GintPrecisionController controller;

    controller.set_mode("mix");
    controller.reset_for_new_scf();
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp32);

    EXPECT_TRUE(controller.update_after_iteration(9.0e-5, 1.0e-7));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp64);
}

TEST(GintPrecisionControllerTest, DefaultModeStartsAndStaysFp64)
{
    GintPrecisionController controller;

    controller.reset_for_new_scf();
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp64);

    EXPECT_FALSE(controller.update_after_iteration(9.0e-5, 1.0e-7));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp64);
}

TEST(GintPrecisionControllerTest, SingleModeStartsAndStaysFp32)
{
    GintPrecisionController controller;

    controller.set_mode("single");
    controller.reset_for_new_scf();
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp32);

    EXPECT_FALSE(controller.update_after_iteration(9.0e-5, 1.0e-7));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp32);
}

TEST(GintPrecisionControllerTest, MixModeLocksFp64AfterSwitch)
{
    GintPrecisionController controller;

    controller.set_mode("mix");
    controller.reset_for_new_scf();
    EXPECT_TRUE(controller.update_after_iteration(9.0e-5, 1.0e-6));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp64);

    // After locking, should return false
    EXPECT_FALSE(controller.update_after_iteration(1.0, 1.0e-6));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp64);
}

TEST(GintPrecisionControllerTest, MixModeReturnsFalseWhenDrhoTooLarge)
{
    GintPrecisionController controller;

    controller.set_mode("mix");
    controller.reset_for_new_scf();
    // drho is large, should not switch yet
    EXPECT_FALSE(controller.update_after_iteration(1.0, 1.0e-7));
    EXPECT_EQ(controller.current_precision(), ModuleGint::GintPrecision::fp32);
}
