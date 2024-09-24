#include "module_base/grid/radial.h"

#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#ifdef __MPI
#include <mpi.h>
#endif

using namespace Grid::Radial;
using Func_t = std::function<double(double)>;

/**
 * This test briefly checks various radial quadrature schemes by comparing
 * their numerical results with analytical values on a few simple functions.
 *
 * The number of grid points and scaling factor are not carefully selected
 * and the test is not exhaustive. It is just a sanity check.
 *
 */

const double pi = std::acos(-1.0);

// tested functions and their analytical integrals
// for \int_0^{+\infty} dr r^2 f(r)
std::vector<std::pair<Func_t, double>> test_func_ref = {
    {
        [](double r) {
            return std::exp(-0.3 * r * r) + std::exp(-3.0 * r * r);
        },
        0.25 * std::sqrt(pi) * (std::pow(0.3, -1.5) + std::pow(3.0, -1.5))
    },
    {
        [](double r) {
            return r * (std::exp(-0.3 * r * r) + std::exp(-3.0 * r * r));
        },
        0.5 / (0.3 * 0.3) + 0.5 / (3.0 * 3.0)
    },
    {
        [](double r) {
            return r * r * (std::exp(-0.3 * r * r) + std::exp(-3.0 * r * r));
        },
        0.375 * std::sqrt(pi) * (std::pow(0.3, -2.5) + std::pow(3.0, -2.5))
    },
};


double quadrature(const Func_t& f, int n, double* r, double* w) {
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        res += w[i] * f(r[i]);
    }
    return res;
}


TEST(RadialTest, Baker) {
    int nbase = 30;
    int mult = 2;
    double R = 2.0;
    std::vector<double> r, w;
    baker(nbase, R, r, w, mult);

    for (auto& t : test_func_ref) {
        double res = quadrature(t.first, r.size(), r.data(), w.data());
        EXPECT_NEAR(res, t.second, 1.0e-5);
    }
}


int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
