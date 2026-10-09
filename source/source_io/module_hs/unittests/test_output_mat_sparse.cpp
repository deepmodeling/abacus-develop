/**
 * @file test_output_mat_sparse.cpp
 * @brief Unit test for output_mat_sparse.h: MatSparseOutputOptions defaults.
 *
 * The struct is the single source of truth for the sparse-matrix output
 * switches; the test pins the legacy default values so an accidental
 * default change is caught at compile-and-run time instead of silently
 * changing output files.
 */
#include <gtest/gtest.h>

#include "source_io/module_hs/output_mat_sparse.h"

TEST(OutputMatSparse, OptionsKeepLegacyDefaults)
{
    ModuleIO::MatSparseOutputOptions options;
    EXPECT_FALSE(options.out_mat_dh);
    EXPECT_FALSE(options.out_mat_ds);
    EXPECT_FALSE(options.out_mat_t);
    EXPECT_FALSE(options.out_mat_r);
    EXPECT_EQ(options.dh_precision, 16);
    EXPECT_EQ(options.ds_precision, 16);
    EXPECT_EQ(options.t_precision, 16);
    EXPECT_EQ(options.r_precision, 16);
    EXPECT_DOUBLE_EQ(options.sparse_threshold, 1e-10);
    EXPECT_FALSE(options.binary);
}
