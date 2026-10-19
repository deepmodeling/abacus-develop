#include "source_base/module_device/device.h"
#include "source_lcao/center2orb_orb11.h"
#include "source_lcao/center2orb_orb21.h"

#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <vector>

namespace
{
Numerical_Orbital_Lm make_orbital(const int l, const int index, const int nr, const int nk)
{
    const double dr = 0.02;
    const double cutoff = (nr - 1) * dr;
    std::vector<double> radial(nr);
    std::vector<double> rab(nr, dr);
    std::vector<double> values(nr);
    for (int ir = 0; ir < nr; ++ir)
    {
        const double r = ir * dr;
        radial[ir] = r;
        values[ir] = std::pow(r, l) * std::pow(cutoff - r, 2) * std::exp(-(0.7 + 0.03 * index) * r);
    }
    Numerical_Orbital_Lm orbital;
    const double* rab_ptr = rab.data();
    const double* radial_ptr = radial.data();
    const double* values_ptr = values.data();
    orbital.set_orbital_info("X",
                             0,
                             l,
                             index,
                             nr,
                             rab_ptr,
                             radial_ptr,
                             Numerical_Orbital_Lm::Psi_Type::Psi,
                             values_ptr,
                             nk,
                             0.1,
                             0.01,
                             false,
                             true,
                             true);
    return orbital;
}

void compare_pairs(const bool use_gpu)
{
    // More than 128 same-grid pairs exercise full and partial batches.
    std::vector<Numerical_Orbital_Lm> basis;
    for (int l = 0; l <= 3; ++l)
    {
        basis.push_back(make_orbital(l, l, 101, 65));
    }
    basis.push_back(make_orbital(2, 4, 81, 65));
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(3);
    gaunt.init_Gaunt(3);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(7, 405, 65);
    std::vector<std::unique_ptr<Center2_Orb::Orb11>> reference;
    std::vector<std::unique_ptr<Center2_Orb::Orb11>> candidate;
    std::vector<Center2_Orb::Orb11*> tasks;
    for (int i = 0; i < 140; ++i)
    {
        const int left = i < 130 ? 1 : i % 5;
        const int right = i < 130 ? 2 : (i + 2) % 5;
        reference.emplace_back(new Center2_Orb::Orb11(basis[left], basis[right], &bessel, gaunt));
        candidate.emplace_back(new Center2_Orb::Orb11(basis[left], basis[right], &bessel, gaunt));
        reference.back()->init_radial_table();
        tasks.push_back(candidate.back().get());
    }
    Center2_Orb::Orb11::init_radial_table(tasks, use_gpu);
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (int i = 0; i < 140; ++i)
    {
        SCOPED_TRACE(i);
        const int left = i < 130 ? 1 : i % 5;
        const int right = i < 130 ? 2 : (i + 2) % 5;
        const int last_left = 2 * basis[left].getL();
        const int last_right = 2 * basis[right].getL();
        for (double distance: {0.0, 0.037, 0.73, 1.91, 3.87, 4.2})
        {
            const double x = 0.8 * distance;
            const double y = 0.36 * distance;
            const double z = 0.48 * distance;
            const ModuleBase::Vector3<double> position(x, y, z);
            for (int ma: {0, last_left})
            {
                for (int mb: {0, last_right})
                {
                    const double expected = reference[i]->cal_overlap(origin, position, ma, mb);
                    const double actual = candidate[i]->cal_overlap(origin, position, ma, mb);
                    const double tolerance = 2e-11 * (1.0 + std::abs(expected));
                    EXPECT_NEAR(actual, expected, tolerance);
                    const auto grad_expected = reference[i]->cal_grad_overlap(origin, position, ma, mb);
                    const auto grad_actual = candidate[i]->cal_grad_overlap(origin, position, ma, mb);
                    const double error = (grad_actual - grad_expected).norm();
                    const double grad_tolerance = 2e-10 * (1.0 + grad_expected.norm());
                    EXPECT_LT(error, grad_tolerance);
                }
            }
        }
    }
}

void compare_products(const bool use_gpu)
{
    std::vector<Numerical_Orbital_Lm> basis;
    for (int l = 0; l <= 3; ++l)
    {
        basis.push_back(make_orbital(l, l, 101, 65));
    }
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(3);
    gaunt.init_Gaunt(3);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(7, 405, 65);
    std::vector<std::unique_ptr<Center2_Orb::Orb21>> reference;
    std::vector<std::unique_ptr<Center2_Orb::Orb21>> candidate;
    std::vector<Center2_Orb::Orb21*> tasks;
    for (int i = 0; i < 4; ++i)
    {
        reference.emplace_back(new Center2_Orb::Orb21(basis[1], basis[2], basis[i], &bessel, gaunt));
        candidate.emplace_back(new Center2_Orb::Orb21(basis[1], basis[2], basis[i], &bessel, gaunt));
        reference.back()->init_radial_table();
        tasks.push_back(candidate.back().get());
    }
    Center2_Orb::Orb21::init_radial_table(tasks, use_gpu);
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (int i = 0; i < 4; ++i)
    {
        for (double distance: {0.0, 0.037, 0.73, 1.91, 3.87})
        {
            const double x = 0.8 * distance;
            const double y = 0.36 * distance;
            const double z = 0.48 * distance;
            const ModuleBase::Vector3<double> position(x, y, z);
            const double expected = reference[i]->cal_overlap(origin, position, 1, 2, i);
            const double actual = candidate[i]->cal_overlap(origin, position, 1, 2, i);
            const double tolerance = 2e-11 * (1.0 + std::abs(expected));
            EXPECT_NEAR(actual, expected, tolerance);
            const auto grad_expected = reference[i]->cal_grad_overlap(origin, position, 1, 2, i);
            const auto grad_actual = candidate[i]->cal_grad_overlap(origin, position, 1, 2, i);
            const double error = (grad_actual - grad_expected).norm();
            const double grad_tolerance = 2e-10 * (1.0 + grad_expected.norm());
            EXPECT_LT(error, grad_tolerance);
        }
    }
}

void compare_high_l_products(const bool use_gpu)
{
    // Heavy-element auxiliary bases include l=6; multiplying by an f orbital
    // requires product angular momenta up to 9 and Bessel functions up to 13.
    const auto left = make_orbital(6, 0, 101, 65);
    const auto middle = make_orbital(3, 1, 101, 65);
    const auto right = make_orbital(3, 2, 81, 65);
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(9);
    gaunt.init_Gaunt(9);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(13, 405, 65);
    Center2_Orb::Orb21 reference(left, middle, right, &bessel, gaunt);
    Center2_Orb::Orb21 candidate(left, middle, right, &bessel, gaunt);
    reference.init_radial_table();
    const std::vector<Center2_Orb::Orb21*> tasks{&candidate};
    Center2_Orb::Orb21::init_radial_table(tasks, use_gpu);
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (double distance: {0.0, 0.037, 0.73, 1.91, 3.57, 3.8})
    {
        const double x = 0.8 * distance;
        const double y = 0.36 * distance;
        const double z = 0.48 * distance;
        const ModuleBase::Vector3<double> position(x, y, z);
        for (int ma: {0, 6, 12})
        {
            for (int mb: {0, 3, 6})
            {
                for (int mc: {0, 3, 6})
                {
                    const double expected = reference.cal_overlap(origin, position, ma, mb, mc);
                    const double actual = candidate.cal_overlap(origin, position, ma, mb, mc);
                    const double tolerance = 2e-11 * (1.0 + std::abs(expected));
                    EXPECT_NEAR(actual, expected, tolerance);
                    const auto grad_expected = reference.cal_grad_overlap(origin, position, ma, mb, mc);
                    const auto grad_actual = candidate.cal_grad_overlap(origin, position, ma, mb, mc);
                    const double error = (grad_actual - grad_expected).norm();
                    const double grad_tolerance = 2e-10 * (1.0 + grad_expected.norm());
                    EXPECT_LT(error, grad_tolerance);
                }
            }
        }
    }
}

void compare_radial_tiles(const bool use_gpu)
{
    // The paired value/derivative Bessel matrix exceeds the 128 MiB workspace
    // limit. Check interpolation across its radial tile boundary and tail.
    const auto orbital = make_orbital(0, 0, 527, 4097);
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(0);
    gaunt.init_Gaunt(0);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(1, 2109, 4097);
    Center2_Orb::Orb11 reference(orbital, orbital, &bessel, gaunt);
    Center2_Orb::Orb11 candidate(orbital, orbital, &bessel, gaunt);
    reference.init_radial_table();
    const std::vector<Center2_Orb::Orb11*> tasks{&candidate};
    Center2_Orb::Orb11::init_radial_table(tasks, use_gpu);
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    for (double distance: {0.0, 0.037, 1.91, 20.465, 20.475, 21.035, 21.1})
    {
        const ModuleBase::Vector3<double> position(distance, 0.0, 0.0);
        const double expected = reference.cal_overlap(origin, position, 0, 0);
        const double actual = candidate.cal_overlap(origin, position, 0, 0);
        const double tolerance = 2e-11 * (1.0 + std::abs(expected));
        EXPECT_NEAR(actual, expected, tolerance);
        const auto grad_expected = reference.cal_grad_overlap(origin, position, 0, 0);
        const auto grad_actual = candidate.cal_grad_overlap(origin, position, 0, 0);
        const double error = (grad_actual - grad_expected).norm();
        const double grad_tolerance = 2e-10 * (1.0 + grad_expected.norm());
        EXPECT_LT(error, grad_tolerance);
    }
}
} // namespace

TEST(Center2Orb11, BatchedCpuMatchesIndividualTables)
{
    compare_pairs(false);
    compare_products(false);
    compare_high_l_products(false);
    compare_radial_tiles(false);
}

TEST(Center2Orb21, ProductGradientMatchesFiniteDifference)
{
    const auto left = make_orbital(2, 0, 101, 65);
    const auto middle = make_orbital(1, 1, 101, 65);
    const auto right = make_orbital(3, 2, 101, 65);
    ORB_gaunt_table gaunt;
    gaunt.init_Gaunt_CH(3);
    gaunt.init_Gaunt(3);
    ModuleBase::Sph_Bessel_Recursive::D2 bessel;
    bessel.set_dx(0.001);
    bessel.cal_jlx(7, 405, 65);
    Center2_Orb::Orb21 orbital(left, middle, right, &bessel, gaunt);
    orbital.init_radial_table();
    const ModuleBase::Vector3<double> origin(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> position(0.813, 0.321, 0.473);
    const double step = 1e-5;
    for (int ma: {0, 2, 4})
    {
        for (int mb: {0, 1, 2})
        {
            for (int mc: {0, 3, 6})
            {
                const auto gradient = orbital.cal_grad_overlap(origin, position, ma, mb, mc);
                auto plus = position;
                auto minus = position;
                plus.x += step;
                minus.x -= step;
                const double upper = orbital.cal_overlap(origin, plus, ma, mb, mc);
                const double lower = orbital.cal_overlap(origin, minus, ma, mb, mc);
                const double finite_difference = (upper - lower) / (2 * step);
                const double tolerance = 2e-7 * (1.0 + std::abs(finite_difference));
                EXPECT_NEAR(gradient.x, finite_difference, tolerance);
            }
        }
    }
}

#ifdef __CUDA
TEST(Center2Orb11, BatchedGpuMatchesIndividualTables)
{
    const bool available = base_device::information::probe_gpu_availability();
    if (!available)
    {
        GTEST_SKIP() << "A CUDA device is required";
    }
    compare_pairs(true);
    compare_products(true);
    compare_high_l_products(true);
    compare_radial_tiles(true);
}
#endif
