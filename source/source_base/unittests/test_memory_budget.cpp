#include "source_base/memory_budget.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <limits>

TEST(MemoryBudget, EstimatorMatchesIssue8342Dimensions)
{
    // NLOCAL ~ 377884 split over 56 ranks gives a local block of ~2.55e9
    // elements; gamma-only genelpa keeps five such dense matrices.
    const std::size_t nrloc = 50500;
    const std::size_t ncloc = 50500;
    const std::size_t local = nrloc * ncloc;
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(nrloc, ncloc, 73820, 1, 1, true, true, true, 8);

    EXPECT_EQ(estimate.hk_bytes, local * 8);
    EXPECT_EQ(estimate.sk_bytes, local * 8);
    EXPECT_EQ(estimate.dmk_bytes, local * 8);
    EXPECT_EQ(estimate.psi_bytes, local * 8);
    EXPECT_EQ(estimate.elpa_work_bytes, local * 8);
    EXPECT_EQ(estimate.per_rank_total_bytes, local * 8 * 5);
    EXPECT_GT(estimate.per_rank_total_bytes * 56, 1500000000000ULL);
}

TEST(MemoryBudget, EstimatorHandlesMultiKAndNonGamma)
{
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(1000, 1000, 100, 4, 1, false, false, false, 8);

    EXPECT_EQ(estimate.hk_bytes, 8 * 1000 * 1000);
    EXPECT_EQ(estimate.sk_bytes, 8 * 1000 * 1000);
    EXPECT_EQ(estimate.dmk_bytes, 4 * 8 * 1000 * 1000);
    EXPECT_EQ(estimate.psi_bytes, 4 * 100 * 1000 * 8);
    EXPECT_EQ(estimate.elpa_work_bytes, 0ULL);
}

TEST(MemoryBudget, EstimatorSaturatesOnOverflow)
{
    const ModuleBase::LcaoDenseEstimate estimate
        = ModuleBase::estimate_lcao_dense(std::numeric_limits<std::size_t>::max(),
                                          std::numeric_limits<std::size_t>::max(),
                                          1,
                                          1,
                                          1,
                                          true,
                                          true,
                                          false,
                                          16);
    EXPECT_EQ(estimate.hk_bytes, std::numeric_limits<std::size_t>::max());
    EXPECT_EQ(estimate.per_rank_total_bytes, std::numeric_limits<std::size_t>::max());
}

TEST(MemoryBudget, GammaDirectSolverSetMatchesSetupPsi)
{
    // Must match the predicate in source_psi/setup_psi.cpp exactly: these
    // solvers keep Psi as a full local N x N block under gamma_only.
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("genelpa"));
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("elpa"));
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("lapack"));
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("pexsi"));
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("cusolver"));
    EXPECT_TRUE(ModuleBase::is_gamma_direct_ks_solver("cusolvermp"));
    // scalapack_gvx / dav / cg keep Psi as N x NBANDS, so they are not
    // gamma-direct even though scalapack_gvx is column-major.
    EXPECT_FALSE(ModuleBase::is_gamma_direct_ks_solver("scalapack_gvx"));
    EXPECT_FALSE(ModuleBase::is_gamma_direct_ks_solver("dav"));
    EXPECT_FALSE(ModuleBase::is_gamma_direct_ks_solver("cg"));
    EXPECT_FALSE(ModuleBase::is_gamma_direct_ks_solver(""));
}

TEST(MemoryBudget, ElpaSolverSet)
{
    EXPECT_TRUE(ModuleBase::is_elpa_ks_solver("genelpa"));
    EXPECT_TRUE(ModuleBase::is_elpa_ks_solver("elpa"));
    EXPECT_FALSE(ModuleBase::is_elpa_ks_solver("lapack"));
    EXPECT_FALSE(ModuleBase::is_elpa_ks_solver("cusolvermp"));
    EXPECT_FALSE(ModuleBase::is_elpa_ks_solver(""));
}

TEST(MemoryBudget, FormatMemoryBytesScalesUnits)
{
    EXPECT_EQ(ModuleBase::format_memory_bytes(0), "0.00 B");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1023), "1023.00 B");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1024), "1.00 KB");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1536), "1.50 KB");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1024ULL * 1024), "1.00 MB");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1024ULL * 1024 * 1024), "1.00 GB");
    EXPECT_EQ(ModuleBase::format_memory_bytes(1024ULL * 1024 * 1024 * 1024), "1.00 TB");
}
