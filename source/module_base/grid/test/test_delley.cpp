#include "module_base/grid/delley.h"
#include "module_base/ylm.h"

#include "gtest/gtest.h"
#include <random>
#ifdef __MPI
#include <mpi.h>
#endif

using namespace Grid::Angular;

// mock the function to prevent unnecessary dependency
namespace ModuleBase {
void WARNING_QUIT(const std::string&, const std::string&) {}
}

class DelleyTest: public ::testing::Test {
protected:
    void randgen(int lmax, std::vector<double>& coef);
    const double tol = 1e-12;
};


void DelleyTest::randgen(int lmax, std::vector<double>& coef) {
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


TEST_F(DelleyTest, NumGrid) {
    int lmax = 5;
    int ngrid = ngrid_delley(lmax);
    EXPECT_EQ(lmax, 17);
    EXPECT_EQ(ngrid, 110);

    lmax = 17;
    ngrid = ngrid_delley(lmax);
    EXPECT_EQ(lmax, 17);
    EXPECT_EQ(ngrid, 110);

    lmax = 20;
    ngrid = ngrid_delley(lmax);
    EXPECT_EQ(lmax, 23);
    EXPECT_EQ(ngrid, 194);

    lmax = 59;
    ngrid = ngrid_delley(lmax);
    EXPECT_EQ(lmax, 59);
    EXPECT_EQ(ngrid, 1202);

    lmax = 60;
    ngrid = ngrid_delley(lmax);
    EXPECT_EQ(lmax, 60);
    EXPECT_EQ(ngrid, -1);
}


TEST_F(DelleyTest, Accuracy) {
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
     * Delley's grid is exact up to floating point errors.
     *
     */
    std::vector<double> grid, weight, coef;

    for (int grid_lmax = 17; grid_lmax < 60; grid_lmax +=6) {
        delley(grid_lmax, grid, weight);
        int func_lmax = grid_lmax / 2;
        randgen(func_lmax, coef);

        double val = 0.0;
        std::vector<double> ylm_real;
        for (size_t i = 0; i < weight.size(); i++) {
            ModuleBase::Ylm::sph_harm(func_lmax,
                    grid[3*i], grid[3*i+1], grid[3*i+2], ylm_real);
            double tmp = 0.0;
            for (size_t j = 0; j < coef.size(); ++j) {
                tmp += coef[j] * ylm_real[j];
            }
            val += weight[i] * tmp * tmp;
        }
        val *= 4.0 * std::acos(-1.0);

        double val_ref = 0.0;
        for (size_t i = 0; i < coef.size(); ++i) {
            val_ref += coef[i] * coef[i];
        }

        double abs_diff = std::abs(val - val_ref);
        EXPECT_LT(abs_diff, tol);
        //printf("order = %2i    val_ref = %8.5f    abs_diff = %8.5e\n",
        //        grid_lmax, val_ref, abs_diff);
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
