#include <gtest/gtest.h>
#include "gint_k.h"

TEST(GintKTest, IntegrationTest) {
    GintK gint;
    gint.initializeParameters();
    double result = gint.computeIntegration();
    EXPECT_NEAR(result, 42.0, 1e-6);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

