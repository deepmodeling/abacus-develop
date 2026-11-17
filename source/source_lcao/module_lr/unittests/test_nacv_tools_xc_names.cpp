#include "../nacv_tools.h"

#include <gtest/gtest.h>

TEST(XcKernelNames, CanonicalizesOnlyActualBuiltinAliases)
{
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("PZ"), "lda");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("SLAPZNOGXNOGC"), "lda");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("SLAPWPBXPBC"), "pbe");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("NOX+NOC"), "none");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("RPA"), "rpa");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("OPT_ORB"), "opt_orb");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("Unknown"), "unknown");
}

TEST(XcKernelNames, PreservesTheExistingCanonicalKernelCatalogue)
{
    EXPECT_TRUE(LR::NacvXC::has_local_xc("lda"));
    EXPECT_TRUE(LR::NacvXC::has_local_xc("pwlda"));
    EXPECT_TRUE(LR::NacvXC::has_local_xc("pbe"));
    EXPECT_EQ(LR::NacvXC::hybrid_xc_list().size(), 9u);
    for (const std::string &name : LR::NacvXC::hybrid_xc_list())
    {
        EXPECT_EQ(LR::NacvXC::has_local_xc(name), name != "hf");
    }
    for (const std::string &name : {"rpa",
             "none",
             "pz",
             "pbesol",
             "scan",
             "blyp_lr",
             "muller"})
    {
        EXPECT_FALSE(LR::NacvXC::has_local_xc(name));
    }
}
