#include "module_base/math_lebedev_laikov.h"
#include "module_base/ylm.h"

#include "gtest/gtest.h"
#include <random>
#ifdef __MPI
#include <mpi.h>
#endif

using ModuleBase::Lebedev_laikov_grid;

// mock the function to prevent unnecessary dependency
namespace ModuleBase {
void WARNING_QUIT(const std::string&, const std::string&) {}
}

class LebedevLaikovTest: public ::testing::Test {
protected:
    void randgen(int lmax, std::vector<double>& coef);
    const double tol = 1e-12;
};


void LebedevLaikovTest::randgen(int lmax, std::vector<double>& coef) {
    coef.resize((lmax + 1) * (lmax + 1));

    // fill coef with uniformly distributed random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (size_t i = 0; i < coef.size(); ++i) {
        coef[i] = dis(gen);
    }

    // normalize the coefficients
    double fac = 0.0;
    for (size_t i = 0; i < coef.size(); ++i) {
        fac += coef[i] * coef[i];
    }

    fac = 1.0 / std::sqrt(fac);
    for (size_t i = 0; i < coef.size(); ++i) {
        coef[i] *= fac;
    }
}


TEST_F(LebedevLaikovTest, Accuracy) {
    /* 
     * Given
     *
     *      f = c[0]*Y00 + c[1]*Y10 + c[2]*Y11 + ...,
     *
     * where c[0], c[1], c[2], ... are some random numbers, the integration
     * of |f|^2 on the unit sphere
     *
     *      \int |f|^2 d\Omega = c[0]^2 + c[1]^2 + c[2]^2 + ... .
     *
     * This test verifies with the above integral that quadrature with
     * Lebedev grid is exact up to floating point errors.
     *
     */

    // (ngrid, lmax)
    std::set<std::pair<int, int>> supported = {
        {6, 3},
        {14, 5},
        {26, 7},
        {38, 9},
        {50, 11},
        {74, 13},
        {86, 15},
        {110, 17},
        {146, 19},
        {170, 21},
        {194, 23},
        {230, 25},
        {266, 27},
        {302, 29},
        {350, 31},
        {434, 35},
        {590, 41},
        {770, 47},
        {974, 53},
        {1202, 59},
        {1454, 65},
        {1730, 71},
        {2030, 77},
        {2354, 83},
        {2702, 89},
        {3074, 95},
        {3470, 101},
        {3890, 107},
        {4334, 113},
        {4802, 119},
        {5294, 125},
        {5810, 131},
    };

    std::vector<double> coef;

    for (auto& grid_info: supported) {
        int ngrid = grid_info.first;
        int grid_lmax = grid_info.second;

        Lebedev_laikov_grid lebgrid(ngrid);
        lebgrid.generate_grid_points();
        
        const double* weight = lebgrid.get_weight();
        const ModuleBase::Vector3<double>* grid = lebgrid.get_grid_coor();

        int func_lmax = grid_lmax / 2;
        randgen(func_lmax, coef);

        double val = 0.0;
        std::vector<double> ylm_real;
        for (int i = 0; i < ngrid; i++) {
            ModuleBase::Ylm::sph_harm(func_lmax,
                    grid[i].x, grid[i].y, grid[i].z, ylm_real);
            double tmp = 0.0;
            for (size_t j = 0; j < coef.size(); ++j) {
                tmp += coef[j] * ylm_real[j];
            }
            val += weight[i] * tmp * tmp;
        }

        double val_ref = 0.0;
        for (size_t i = 0; i < coef.size(); ++i) {
            val_ref += coef[i] * coef[i];
        }

        double abs_diff = std::abs(val - val_ref);
        EXPECT_LT(abs_diff, tol);
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
