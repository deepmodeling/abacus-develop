#include "source_base/module_device/device.h"
#include "source_lcao/center2orb_orb11.h"

#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <vector>

namespace
{
Numerical_Orbital_Lm make_orbital(const int l, const int index, const int nr)
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
                             65,
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
        basis.push_back(make_orbital(l, l, 101));
    }
    basis.push_back(make_orbital(2, 4, 81));
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
} // namespace

TEST(Center2Orb11, BatchedCpuMatchesIndividualTables)
{
    compare_pairs(false);
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
}
#endif
