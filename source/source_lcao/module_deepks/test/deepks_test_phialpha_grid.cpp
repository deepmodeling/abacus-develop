#include "deepks_test_runner.h"

#include "source_lcao/module_deepks/deepks_iterate.h"
#include "source_lcao/module_rt/snap_phialpha_half_tddft.h"
#include "source_lcao/module_rt/snap_projector_half_tddft.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <type_traits>

template <typename T>
void test_deepks<T>::check_phialpha_grid_zero_field()
{
    struct ComparisonStats
    {
        double max_real_diff = 0.0;
        double max_imag_abs = 0.0;
        double max_reference_abs = 0.0;
        double max_relative_diff = 0.0;
        double reference_at_max_diff = 0.0;
        int compared = 0;
        bool shape_mismatch = false;
    };

    const ModuleBase::Vector3<double> zero_A(0.0, 0.0, 0.0);
    const auto compare_grid = [&](const int radial_grid_num, const int lebedev_grid_points) {
        ComparisonStats stats;
        module_rt::SnapIntegrationOptions options;
        options.radial_grid_num = radial_grid_num;
        options.lebedev_grid_points = lebedev_grid_points;

        DeePKS_domain::iterate_ad1(
            ucell,
            Test_Deepks::GridD,
            ORB,
            false,
            [&](const int iat,
                const ModuleBase::Vector3<double>& tau0,
                const int ibt,
                const ModuleBase::Vector3<double>& tau1,
                const int start,
                const int nw_tot,
                ModuleBase::Vector3<int> dR) {
                const int T1 = ucell.iat2it[ibt];
                const Atom* atom1 = &ucell.atoms[T1];

                auto all_indexes = ParaO.get_indexes_row(ibt);
                auto col_indexes = ParaO.get_indexes_col(ibt);
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());

                for (size_t iw1l = 0; iw1l < all_indexes.size(); iw1l += this->npol)
                {
                    const int iw1 = all_indexes[iw1l] / this->npol;
                    const int L1 = atom1->iw2l[iw1];
                    const int N1 = atom1->iw2n[iw1];
                    const int m1 = atom1->iw2m[iw1];
                    const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                    std::vector<std::vector<std::complex<double>>> grid_nlm;
                    module_rt::snap_phialpha_half_tddft(ORB,
                                                        grid_nlm,
                                                        tau1 * ucell.lat0,
                                                        T1,
                                                        L1,
                                                        m1,
                                                        N1,
                                                        tau0 * ucell.lat0,
                                                        zero_A,
                                                        false,
                                                        options);

                    std::vector<std::vector<double>> tci_nlm;
                    const int T0_fixed = 0;
                    overlap_orb_alpha_.snap(T1,
                                            L1,
                                            N1,
                                            M1,
                                            T0_fixed,
                                            (tau0 - tau1) * ucell.lat0,
                                            false,
                                            tci_nlm);

                    if (grid_nlm.empty() || tci_nlm.empty() || grid_nlm[0].size() != tci_nlm[0].size())
                    {
                        stats.shape_mismatch = true;
                        return;
                    }

                    for (size_t i = 0; i < grid_nlm[0].size(); ++i)
                    {
                        const double reference_abs = std::abs(tci_nlm[0][i]);
                        const double real_diff = std::abs(grid_nlm[0][i].real() - tci_nlm[0][i]);
                        if (real_diff > stats.max_real_diff)
                        {
                            stats.max_real_diff = real_diff;
                            stats.reference_at_max_diff = tci_nlm[0][i];
                        }
                        stats.max_imag_abs = std::max(stats.max_imag_abs, std::abs(grid_nlm[0][i].imag()));
                        stats.max_reference_abs = std::max(stats.max_reference_abs, reference_abs);
                        if (reference_abs > 1.0e-8)
                        {
                            stats.max_relative_diff = std::max(stats.max_relative_diff, real_diff / reference_abs);
                        }
                        ++stats.compared;
                    }
                }
            });

        const char* instance = std::is_same<T, double>::value ? "gamma" : "multik";
        std::cout << std::scientific << std::setprecision(12) << "phialpha " << instance << " grid "
                  << radial_grid_num << "x" << lebedev_grid_points
                  << ": max abs error = " << stats.max_real_diff
                  << " (reference = " << stats.reference_at_max_diff << ")"
                  << ", max reference = " << stats.max_reference_abs
                  << ", max imaginary magnitude = " << stats.max_imag_abs
                  << ", max relative error (|reference| > 1e-8) = " << stats.max_relative_diff
                  << std::defaultfloat << std::endl;

        EXPECT_FALSE(stats.shape_mismatch) << "phialpha grid and two-center integration output shapes differ";
        EXPECT_GT(stats.compared, 0) << "No phialpha grid entries were compared";
        EXPECT_LE(stats.max_imag_abs, 1.0e-14) << "max reference abs = " << stats.max_reference_abs;
        return stats;
    };

    const ComparisonStats default_grid = compare_grid(140, 110);
    const ComparisonStats dense_radial_grid = compare_grid(280, 110);
    const ComparisonStats dense_angular_grid = compare_grid(140, 590);

    const bool is_gamma = std::is_same<T, double>::value;
    const double default_grid_tolerance = is_gamma ? 6.0e-5 : 1.0e-5;
    const double dense_angular_grid_tolerance = is_gamma ? 5.0e-6 : 4.0e-6;

    // The 110-point angular rule limits both radial-grid cases. The 590-point
    // rule exposes the lower error reached by the corrected interpolation.
    EXPECT_LE(default_grid.max_real_diff, default_grid_tolerance);
    EXPECT_LE(dense_radial_grid.max_real_diff, default_grid_tolerance);
    EXPECT_NEAR(dense_radial_grid.max_real_diff, default_grid.max_real_diff, 2.0e-10);
    EXPECT_LE(dense_angular_grid.max_real_diff, dense_angular_grid_tolerance);
    EXPECT_LE(dense_angular_grid.max_real_diff, 0.5 * default_grid.max_real_diff);
}

template void test_deepks<double>::check_phialpha_grid_zero_field();
template void test_deepks<std::complex<double>>::check_phialpha_grid_zero_field();

template <typename T>
void run_deepks_unit_phialpha_grid_zero_field(test_deepks<T>& test)
{
    test.check_phialpha_grid_zero_field();
}

template void run_deepks_unit_phialpha_grid_zero_field<double>(test_deepks<double>& test);
template void run_deepks_unit_phialpha_grid_zero_field<std::complex<double>>(test_deepks<std::complex<double>>& test);
