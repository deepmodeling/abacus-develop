#include "../libxc_abacus.h"

#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#ifdef __LIBXC

TEST(LibxcSanitizer, FullPolarizationUsesTheWeightedEnergyDerivative)
{
    const std::array<int, 2> functional_ids = {{XC_LDA_X, XC_LDA_C_PZ}};
    const double density = 0.45;

    for (std::size_t ifunc = 0; ifunc < functional_ids.size(); ++ifunc)
    {
        xc_func_type func;
        ASSERT_EQ(xc_func_init(&func, functional_ids[ifunc], XC_POLARIZED), 0);
        xc_func_set_dens_threshold(&func, 1.0e-6);

        const auto evaluate = [&func](const double rho_up,
                                      const double rho_down,
                                      XC_Functional_Libxc::LibxcWeightedDerivatives* const weighted) {
            const std::vector<double> rho = {rho_up, rho_down};
            const std::vector<double> mask = {1.0, 1.0};
            std::vector<double> exc(1, 0.0);
            std::vector<double> vrho(2, 0.0);
            xc_lda_exc_vxc(&func, 1, rho.data(), exc.data(), vrho.data());
            if (weighted != nullptr)
            {
                *weighted = XC_Functional_Libxc::make_libxc_weighted_derivatives(func,
                                                                                 2,
                                                                                 1,
                                                                                 mask,
                                                                                 rho,
                                                                                 std::vector<double>(),
                                                                                 exc,
                                                                                 vrho,
                                                                                 std::vector<double>());
            }
            return (rho_up + rho_down) * exc[0];
        };

        XC_Functional_Libxc::LibxcWeightedDerivatives weighted;
        const double energy = evaluate(density, 0.0, &weighted);
        ASSERT_EQ(weighted.drho.size(), 2U);
        EXPECT_DOUBLE_EQ(weighted.energy_sum, energy);

        const double steps[] = {1.0e-3, 5.0e-4, 2.5e-4, 1.25e-4};
        std::array<double, 4> errors = {{0.0, 0.0, 0.0, 0.0}};
        for (std::size_t ieps = 0; ieps < 4; ++ieps)
        {
            const double step = steps[ieps];
            const double finite_difference
                = (evaluate(density + step, 0.0, nullptr) - evaluate(density - step, 0.0, nullptr)) / (2.0 * step);
            errors[ieps] = std::abs(weighted.drho[0] - finite_difference);
            EXPECT_LE(errors[ieps], 2.0e-7 * std::max(1.0, std::abs(weighted.drho[0])))
                << "functional_id=" << functional_ids[ifunc] << ", step=" << step;
        }
        EXPECT_LE(errors[1], 0.4 * errors[0] + 1.0e-12);
        EXPECT_LE(errors[2], 0.4 * errors[1] + 1.0e-12);

        const double inactive_density = 0.5 * func.dens_threshold;
        const double inactive_step = 0.2 * func.dens_threshold;
        XC_Functional_Libxc::LibxcWeightedDerivatives inactive_weighted;
        evaluate(density, inactive_density, &inactive_weighted);
        const double inactive_finite_difference = (evaluate(density, inactive_density + inactive_step, nullptr)
                                                   - evaluate(density, inactive_density - inactive_step, nullptr))
                                                  / (2.0 * inactive_step);
        EXPECT_NEAR(inactive_finite_difference,
                    inactive_weighted.drho[1],
                    2.0e-8 * std::max(1.0, std::abs(inactive_finite_difference)));
        xc_func_end(&func);
    }
}

TEST(LibxcSanitizer, GgaSigmaReverseMatchesTheWeightedEnergy)
{
    xc_func_type func;
    ASSERT_EQ(xc_func_init(&func, XC_GGA_C_PBE, XC_POLARIZED), 0);
    xc_func_set_dens_threshold(&func, 1.0e-6);
    xc_func_set_sigma_threshold(&func, 1.0e-2);

    const std::vector<double> mask = {1.0, 0.0};
    const std::vector<double> density = {0.40, 0.20};
    const double sigma_floor = func.sigma_threshold * func.sigma_threshold;
    const std::array<std::array<double, 3>, 5> sigma_states = {{{{0.040, 0.010, 0.030}},
                                                                {{0.5 * sigma_floor, 0.0, 0.030}},
                                                                {{0.040, 0.200, 0.030}},
                                                                {{0.040, -0.200, 0.030}},
                                                                {{0.5 * sigma_floor, 0.200, 0.030}}}};

    const auto evaluate = [&func, &mask](const std::vector<double>& rho,
                                         const std::vector<double>& sigma,
                                         XC_Functional_Libxc::LibxcWeightedDerivatives* const weighted) {
        std::vector<double> exc(1, 0.0);
        std::vector<double> vrho(2, 0.0);
        std::vector<double> vsigma(3, 0.0);
        xc_gga_exc_vxc(&func, 1, rho.data(), sigma.data(), exc.data(), vrho.data(), vsigma.data());
        if (weighted != nullptr)
        {
            *weighted
                = XC_Functional_Libxc::make_libxc_weighted_derivatives(func, 2, 1, mask, rho, sigma, exc, vrho, vsigma);
        }
        return (mask[0] * rho[0] + mask[1] * rho[1]) * exc[0];
    };

    for (std::size_t icase = 0; icase < sigma_states.size(); ++icase)
    {
        std::vector<double> sigma(sigma_states[icase].begin(), sigma_states[icase].end());
        XC_Functional_Libxc::LibxcWeightedDerivatives weighted;
        const double energy = evaluate(density, sigma, &weighted);
        EXPECT_DOUBLE_EQ(weighted.energy_sum, energy);
        ASSERT_EQ(weighted.drho.size(), 2U);
        ASSERT_EQ(weighted.dsigma.size(), 3U);

        if (icase == 0)
        {
            for (int component = 0; component < 2; ++component)
            {
                std::vector<double> perturbed_density = density;
                const double step = 1.0e-6;
                perturbed_density[component] += step;
                const double energy_plus = evaluate(perturbed_density, sigma, nullptr);
                perturbed_density[component] -= 2.0 * step;
                const double energy_minus = evaluate(perturbed_density, sigma, nullptr);
                const double finite_difference = (energy_plus - energy_minus) / (2.0 * step);
                EXPECT_NEAR(finite_difference,
                            weighted.drho[component],
                            2.0e-7 * std::max(1.0, std::abs(weighted.drho[component])));
            }
        }

        if (icase == 1 || icase == 4)
        {
            EXPECT_DOUBLE_EQ(weighted.dsigma[0], 0.0);
        }
        if (icase >= 2)
        {
            EXPECT_DOUBLE_EQ(weighted.dsigma[1], 0.0);
        }

        for (int component = 0; component < 3; ++component)
        {
            const double step = 1.0e-6;
            sigma[component] += step;
            const double energy_plus = evaluate(density, sigma, nullptr);
            sigma[component] -= 2.0 * step;
            const double energy_minus = evaluate(density, sigma, nullptr);
            sigma[component] += step;
            const double finite_difference = (energy_plus - energy_minus) / (2.0 * step);
            EXPECT_NEAR(finite_difference,
                        weighted.dsigma[component],
                        2.0e-7 * std::max(1.0, std::abs(weighted.dsigma[component])))
                << "case=" << icase << ", sigma component=" << component;
        }
    }
    xc_func_end(&func);
}

TEST(LibxcSanitizer, UnpolarizedSelfSigmaReverseMatchesTheWeightedEnergy)
{
    xc_func_type func;
    ASSERT_EQ(xc_func_init(&func, XC_GGA_C_PBE, XC_UNPOLARIZED), 0);
    xc_func_set_dens_threshold(&func, 1.0e-6);
    xc_func_set_sigma_threshold(&func, 1.0e-2);

    const std::vector<double> mask = {1.0};
    const std::vector<double> density = {0.40};
    const double sigma_floor = func.sigma_threshold * func.sigma_threshold;

    const auto evaluate = [&func, &mask, &density](const double sigma_value,
                                                   XC_Functional_Libxc::LibxcWeightedDerivatives* const weighted) {
        const std::vector<double> sigma = {sigma_value};
        std::vector<double> exc(1, 0.0);
        std::vector<double> vrho(1, 0.0);
        std::vector<double> vsigma(1, 0.0);
        xc_gga_exc_vxc(&func, 1, density.data(), sigma.data(), exc.data(), vrho.data(), vsigma.data());
        if (weighted != nullptr)
        {
            *weighted = XC_Functional_Libxc::make_libxc_weighted_derivatives(func,
                                                                             1,
                                                                             1,
                                                                             mask,
                                                                             density,
                                                                             sigma,
                                                                             exc,
                                                                             vrho,
                                                                             vsigma);
        }
        return mask[0] * density[0] * exc[0];
    };

    const double below_floor_sigma = 0.5 * sigma_floor;
    const double below_floor_step = 0.2 * sigma_floor;
    XC_Functional_Libxc::LibxcWeightedDerivatives below_floor_weighted;
    const double below_floor_energy = evaluate(below_floor_sigma, &below_floor_weighted);
    ASSERT_EQ(below_floor_weighted.dsigma.size(), 1U);
    EXPECT_DOUBLE_EQ(below_floor_weighted.energy_sum, below_floor_energy);
    EXPECT_DOUBLE_EQ(below_floor_weighted.dsigma[0], 0.0);
    const double below_floor_finite_difference = (evaluate(below_floor_sigma + below_floor_step, nullptr)
                                                  - evaluate(below_floor_sigma - below_floor_step, nullptr))
                                                 / (2.0 * below_floor_step);
    EXPECT_NEAR(below_floor_finite_difference, below_floor_weighted.dsigma[0], 1.0e-12);

    const double above_floor_sigma = 0.040;
    const double above_floor_step = 1.0e-6;
    XC_Functional_Libxc::LibxcWeightedDerivatives above_floor_weighted;
    const double above_floor_energy = evaluate(above_floor_sigma, &above_floor_weighted);
    ASSERT_EQ(above_floor_weighted.dsigma.size(), 1U);
    EXPECT_DOUBLE_EQ(above_floor_weighted.energy_sum, above_floor_energy);
    EXPECT_GT(std::abs(above_floor_weighted.dsigma[0]), 1.0e-12);
    const double above_floor_finite_difference = (evaluate(above_floor_sigma + above_floor_step, nullptr)
                                                  - evaluate(above_floor_sigma - above_floor_step, nullptr))
                                                 / (2.0 * above_floor_step);
    EXPECT_NEAR(above_floor_finite_difference,
                above_floor_weighted.dsigma[0],
                2.0e-7 * std::max(1.0, std::abs(above_floor_weighted.dsigma[0])));

    xc_func_end(&func);
}

#endif
