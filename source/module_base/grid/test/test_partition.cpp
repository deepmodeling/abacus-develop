#include "module_base/grid/partition.h"
#include "module_base/grid/radial.h"
#include "module_base/grid/delley.h"
#include "module_base/constants.h"

#include "gtest/gtest.h"
#include <cmath>
#include <array>
#include <numeric>
#include <chrono>
#include <random>
#include <algorithm>

#ifdef __MPI
#include <mpi.h>
#endif

using ModuleBase::PI;
using Vec3 = std::array<double, 3>;

using iclock = std::chrono::high_resolution_clock;
iclock::time_point start;
std::chrono::duration<double> dur;

double norm(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

Vec3 operator-(const Vec3& v1, const Vec3& v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

Vec3 operator+(const Vec3& v1, const Vec3& v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

// |r|^n * exp(-a*|r|^2)
double func_core(const Vec3& r, double a, double n) {
    double rabs = norm(r);
    return std::pow(rabs, n) * std::exp(-a * rabs * rabs);
}

// func_core integrated over all space
double ref_core(double a, double n) {
    double p = 0.5 * (n + 3);
    return 2.0 * PI * std::pow(a, -p) * std::tgamma(p);
}

// the test function is a combination of several func_core
double func(
    const Vec3& r,
    const std::vector<Vec3>& R,
    const std::vector<double>& a,
    const std::vector<double>& n
) {
    double val = 0.0;
    for (size_t i = 0; i < R.size(); i++) {
        val += func_core(r - R[i], a[i], n[i]);
    }
    return val;
}

double ref(const std::vector<double>& a, const std::vector<double>& n) {
    double val = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        val += ref_core(a[i], n[i]);
    }
    return val;
}

// A Param object specifies a test function
struct Param {
    std::vector<Vec3> R;
    std::vector<double> a;
    std::vector<double> n;
};

std::vector<Param> test_params = {
    {
        {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 2.0},
        },
        {0.5, 2.0},
        {0, 0}
    },
    {
        {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 2.0},
            {0.0, 3.0, 0.0},
        },
        {0.5, 2.0, 1.5},
        {1, 2, 0.5}
    },
    {
        {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 3.0},
            {0.0, 3.0, 0.0},
            {9.0, 0.0, 0.0},
        },
        {1.0, 2.0, 1.5, 2.0},
        {2.5, 2, 0.5, 1}
    },
    {
        {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 3.0},
            {0.0, 3.0, 0.0},
            {9.0, 0.0, 0.0},
            {1.0, 1.0, 1.0},
            {2.0, 2.0, 2.0},
            {3.0, 3.0, 3.0},
            {4.0, 4.0, 4.0},
            {5.0, 5.0, 5.0},
            {6.0, 6.0, 6.0},
        },
        {1.0, 2.0, 1.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0},
        {2.5, 2.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
    },
};

std::vector<double> dist_R_R(const std::vector<Vec3>& R) {
    // tabulate dRR[I,J] = || R[I] - R[J] ||
    size_t nR = R.size();
    std::vector<double> dRR(nR*nR, 0.0);
    for (size_t I = 0; I < nR; I++) {
        for (size_t J = I + 1; J < nR; J++) {
            double d = norm(R[I] - R[J]);
            dRR[I*nR + J] = d;
            dRR[J*nR + I] = d;
        }
    }
    return dRR;
}

class PartitionTest: public ::testing::Test {
protected:
    PartitionTest();

    // grid & weight for one-center integration
    std::vector<double> r;
    std::vector<double> w;

    const double tol = 1e-5;
};

PartitionTest::PartitionTest() {
    // angular grid & weight
    std::vector<double> r_ang, w_ang;
    int lmax = 25;
    Grid::Angular::delley(lmax, r_ang, w_ang);

    // radial grid & weight
    std::vector<double> r_rad, w_rad;
    int nrad = 60;
    int Rcut = 7.0;
    int mult = 2;
    Grid::Radial::baker(nrad, Rcut, r_rad, w_rad, mult);

    // complete grid & weight for one-center integration
    size_t ngrid = w_rad.size() * w_ang.size();
    r.resize(3*ngrid);
    w.resize(ngrid);

    size_t ir = 0;
    for (size_t i = 0; i < w_rad.size(); i++) {
        for (size_t j = 0; j < w_ang.size(); j++) {
            r[3*ir] = r_rad[i] * r_ang[3*j];
            r[3*ir+1] = r_rad[i] * r_ang[3*j+1];
            r[3*ir+2] = r_rad[i] * r_ang[3*j+2];
            w[ir] = w_rad[i] * w_ang[j] * 4.0 * PI;
            ++ir;
        }
    }
}


TEST_F(PartitionTest, Becke) {
    dur = dur.zero();
    for (const Param& param : test_params) {
        double val = 0.0;
        double val_ref = ref(param.a, param.n);

        // tabulate || R[I] - R[J] ||
        std::vector<double> dRR(dist_R_R(param.R));

        // all centers are involved
        size_t nR = param.R.size();
        std::vector<int> iR(nR);
        std::iota(iR.begin(), iR.end(), 0);

        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(iR.begin(), iR.end(), g);

        for (size_t I = 0; I < nR; ++I) { // for each center
            for (size_t i = 0; i < w.size(); i++) {
                Vec3 ri = Vec3{r[3*i], r[3*i+1], r[3*i+2]} + param.R[I];

                // tabulate || r - R[J] ||
                std::vector<double> drR(nR);
                for (size_t J = 0; J < nR; ++J) {
                    drR[J] = norm(ri - param.R[J]);
                }

                // partition weight for this grid point
                start = iclock::now();
                double w_part = Grid::Partition::w_becke(
                    drR.size(), drR.data(), dRR.data(),
                    iR.size(), iR.data(), I
                );
                dur += iclock::now() - start;

                val += w_part * w[i] * func(ri, param.R, param.a, param.n);
            }
        }

        EXPECT_NEAR(val, val_ref, tol);
    }
    printf("time elapsed = %8.3e seconds\n", dur.count());
}


TEST_F(PartitionTest, Stratmann) {
    dur = dur.zero();

    for (const Param& param : test_params) {
        double val = 0.0;
        double val_ref = ref(param.a, param.n);

        // tabulate || R[I] - R[J] ||
        std::vector<double> dRR(dist_R_R(param.R));

        // all centers are involved
        size_t nR = param.R.size();
        std::vector<int> iR(nR);
        std::iota(iR.begin(), iR.end(), 0);

        // radii of exclusive zone
        std::vector<double> drR_thr(nR);
        for (size_t I = 0; I < nR; ++I) {
            double dRRmin = 1e100;
            for (size_t J = 0; J < nR; ++J) {
                if (J != I) {
                    dRRmin = std::min(dRRmin, dRR[I*nR + J]);
                }
            }
            drR_thr[I] = 0.5 * (1.0 - Grid::Partition::stratmann_a) * dRRmin;
        }

        for (size_t I = 0; I < nR; ++I) { // for each center
            for (size_t i = 0; i < w.size(); i++) {
                Vec3 ri = Vec3{r[3*i], r[3*i+1], r[3*i+2]} + param.R[I];

                // tabulate || r - R[J] ||
                std::vector<double> drR(nR);
                for (size_t J = 0; J < nR; ++J) {
                    drR[J] = norm(ri - param.R[J]);
                }

                // partition weight for this grid point
                start = iclock::now();
                double w_part = Grid::Partition::w_stratmann(
                    drR.size(), drR.data(), dRR.data(), drR_thr.data(), 
                    iR.size(), iR.data(), I
                );
                dur += iclock::now() - start;

                val += w_part * w[i] * func(ri, param.R, param.a, param.n);
            }
        }

        EXPECT_NEAR(val, val_ref, tol);
    }
    printf("time elapsed = %8.3e seconds\n", dur.count());
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
