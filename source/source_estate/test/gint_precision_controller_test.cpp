#include "gtest/gtest.h"

#include "../module_charge/gint_precision_controller.h"

#include <cstdlib>

namespace
{
class ScopedEnvVar
{
  public:
    ScopedEnvVar(const char* name, const char* value)
        : name_(name)
    {
        const char* old = std::getenv(name_);
        if (old != nullptr)
        {
            had_old_value_ = true;
            old_value_ = old;
        }

        if (value == nullptr)
        {
            unsetenv(name_);
        }
        else
        {
            setenv(name_, value, 1);
        }
    }

    ~ScopedEnvVar()
    {
        if (had_old_value_)
        {
            setenv(name_, old_value_.c_str(), 1);
        }
        else
        {
            unsetenv(name_);
        }
    }

  private:
    const char* name_;
    bool had_old_value_ = false;
    std::string old_value_;
};
}

TEST(GintPrecisionControllerTest, AutoModeSwitchesToFp64AfterTwoQualifiedIterations)
{
    ScopedEnvVar env_guard("ABACUS_GINT_FORCE_CPU_REAL", nullptr);
    GintPrecisionController controller;

    controller.reset_for_new_scf();
    EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);

    controller.update_after_iteration(3, 9.0e-5, 1.0e-6, false, false);
    EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);

    controller.update_after_iteration(4, 8.0e-5, 1.0e-6, false, false);
    EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);
}

TEST(GintPrecisionControllerTest, RuntimeOverrideCanForceFp32OrFp64)
{
    {
        ScopedEnvVar env_guard("ABACUS_GINT_FORCE_CPU_REAL", "fp64");
        GintPrecisionController controller;
        controller.reset_for_new_scf();
        EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);
        controller.update_after_iteration(4, 9.0e-5, 1.0e-6, false, false);
        EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp64);
    }

    {
        ScopedEnvVar env_guard("ABACUS_GINT_FORCE_CPU_REAL", "fp32");
        GintPrecisionController controller;
        controller.reset_for_new_scf();
        EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);
        controller.update_after_iteration(3, 1.0e-4, 1.0e-6, false, false);
        controller.update_after_iteration(4, 9.0e-5, 1.0e-6, false, false);
        EXPECT_EQ(controller.current_config().cpu_internal_real, ModuleGint::GintRealPrecision::fp32);
    }
}
