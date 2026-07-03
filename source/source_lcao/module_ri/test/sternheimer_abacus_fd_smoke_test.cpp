#include "source_lcao/module_ri/sternheimer_abacus_fd_smoke.h"

#include <gtest/gtest.h>
#include <string>

TEST(SternheimerABACUSFDSmoke, FormatsZeroOrderComparisonReport)
{
    ModuleRI::SternheimerABACUSDFTZeroOrderResult result;
    result.grid_data.grid.nx = 2;
    result.grid_data.grid.ny = 3;
    result.grid_data.grid.nz = 4;
    result.grid_data.grid.hx = 0.5;
    result.grid_data.grid.hy = 0.5;
    result.grid_data.grid.hz = 0.5;
    result.grid_data.volume_element = 0.125;
    result.fd_states.eigenvalues = {-0.4, 0.2};
    result.fd_states.residual_norms = {1.0e-10, 2.0e-9};
    result.dft_eigenvalues = {-0.41, 0.25};
    result.dft_occupations = {1.0, 0.0};
    result.hamiltonian_mode = "full_nc_nonlocal";
    result.lanczos_max_subspace_size = 640;

    result.comparison.max_abs_fd_minus_dft = 0.05;
    result.comparison.all_eigenvalues_within_tolerance = false;
    result.comparison.bands = {
        ModuleRI::SternheimerDFTZeroOrderBandComparison{0, -0.4, -0.41, 0.01, 1.0e-10, true},
        ModuleRI::SternheimerDFTZeroOrderBandComparison{1, 0.2, 0.25, -0.05, 2.0e-9, false},
    };

    const std::string report = ModuleRI::format_sternheimer_fd_zero_order_report(result);

    EXPECT_NE(report.find("# ABACUS Sternheimer FD zero-order smoke test"), std::string::npos);
    EXPECT_NE(report.find("grid 2 3 4 size 24 dV 0.125"), std::string::npos);
    EXPECT_NE(report.find("hamiltonian_mode full_nc_nonlocal"), std::string::npos);
    EXPECT_NE(report.find("lanczos_max_subspace_size 640"), std::string::npos);
    EXPECT_NE(report.find("band fd_eigenvalue_Ry dft_eigenvalue_Ry fd_minus_dft_Ry occupation fd_residual_norm eigenvalue_within_tolerance"), std::string::npos);
    EXPECT_NE(report.find("0 -0.4 -0.41 0.01 1 1e-10 yes"), std::string::npos);
    EXPECT_NE(report.find("1 0.2 0.25 -0.05 0 2e-09 no"), std::string::npos);
    EXPECT_NE(report.find("max_abs_fd_minus_dft_Ry 0.05"), std::string::npos);
    EXPECT_NE(report.find("all_eigenvalues_within_tolerance no"), std::string::npos);
}
