#include "../potentials/xc_kxc.h"
#include <xc.h>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

TEST(XcKxc, NoGradientDoesNotRequireKxc)
{
    xc_func_info_type info = {};
    info.name = "missing KXC";
    xc_func_type functional = {};
    functional.info = &info;
    EXPECT_NO_THROW(LR::require_xc_kxc(functional, false));
}

TEST(XcKxc, MissingKxcNamesFunctionalAndBuildOption)
{
    xc_func_info_type info = {};
    info.name = "missing KXC";
    xc_func_type functional = {};
    functional.info = &info;
    try
    {
        LR::require_xc_kxc(functional, true);
        FAIL() << "missing KXC was accepted";
    }
    catch (const std::runtime_error& error)
    {
        const std::string message = error.what();
        const std::string version = xc_version_string();
        EXPECT_NE(message.find(info.name), std::string::npos);
        EXPECT_NE(message.find("-DDISABLE_KXC=OFF"), std::string::npos);
        EXPECT_NE(message.find(version), std::string::npos);
    }
}

TEST(XcKxc, MixedFunctionalChecksComponentEvenWhenParentClaimsKxc)
{
    xc_func_info_type parent_info = {};
    parent_info.flags = XC_FLAGS_HAVE_KXC;
    parent_info.name = "hybrid";
    xc_func_info_type child_info = {};
    child_info.name = "WPBEH component";
    xc_func_type child = {};
    child.info = &child_info;
    xc_func_type* children[] = {&child};
    double coefficients[] = {1.0};
    xc_func_type parent = {};
    parent.info = &parent_info;
    parent.n_func_aux = 1;
    parent.func_aux = children;
    parent.mix_coef = coefficients;
    EXPECT_THROW(LR::require_xc_kxc(parent, true), std::runtime_error);
    child_info.flags = XC_FLAGS_HAVE_KXC;
    EXPECT_NO_THROW(LR::require_xc_kxc(parent, true));
}

TEST(XcKxc, InternalAuxiliaryIsNotAnIndependentMixedContribution)
{
    xc_func_info_type info = {};
    info.flags = XC_FLAGS_HAVE_KXC;
    xc_func_type internal = {};
    xc_func_type* children[] = {&internal};
    xc_func_type functional = {};
    functional.info = &info;
    functional.n_func_aux = 1;
    functional.func_aux = children;
    EXPECT_NO_THROW(LR::require_xc_kxc(functional, true));
}

TEST(XcKxc, RuntimeLdaFlagsDetermineAvailability)
{
    xc_func_type functional = {};
    ASSERT_EQ(xc_func_init(&functional, XC_LDA_X, XC_UNPOLARIZED), 0);
    const xc_func_info_type* info = xc_func_get_info(&functional);
    const int flags = xc_func_info_get_flags(info);
    if ((flags & XC_FLAGS_HAVE_KXC) != 0)
    {
        EXPECT_NO_THROW(LR::require_xc_kxc(functional, true));
    }
    else { EXPECT_THROW(LR::require_xc_kxc(functional, true), std::runtime_error); }
    xc_func_end(&functional);
}
