// g++ test_delley.cpp ../delley.cpp ../../ylm.cpp ../../tool_quit.cpp ../../global_variable.cpp ../../timer.cpp ../../global_file.cpp ../../global_function.cpp -I../../.. -lgtest -o test_delley
#include "module_base/grid/delley.h"

#include "module_base/ylm.h"

#include "gtest/gtest.h"
#include <random>
#ifdef __MPI
#include <mpi.h>
#endif

class DelleyTest: public ::testing::Test {

protected:
    void randgen(int lmax, std::vector<double>& coef);
    const double rel_tol = 1e-12;
    const double abs_tol = 1e-12;
};

void DelleyTest::randgen(int lmax, std::vector<double>& coef) {
    coef.resize(lmax + 1);

    // fill coef with uniformly distributed random numbers
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i = 0; i <= lmax; i++) {
        coef[i] = dis(gen) / std::sqrt(lmax);
    }
}


TEST_F(DelleyTest, NumGrid) {
    int lmax = 5;
    int ngrid = Delley::ngrid(lmax);
    EXPECT_EQ(lmax, 17);
    EXPECT_EQ(ngrid, 110);

    lmax = 17;
    ngrid = Delley::ngrid(lmax);
    EXPECT_EQ(lmax, 17);
    EXPECT_EQ(ngrid, 110);

    lmax = 20;
    ngrid = Delley::ngrid(lmax);
    EXPECT_EQ(lmax, 23);
    EXPECT_EQ(ngrid, 194);

    lmax = 59;
    ngrid = Delley::ngrid(lmax);
    EXPECT_EQ(lmax, 59);
    EXPECT_EQ(ngrid, 1202);

    lmax = 60;
    ngrid = Delley::ngrid(lmax);
    EXPECT_EQ(lmax, 60);
    EXPECT_EQ(ngrid, -1);
}


TEST_F(DelleyTest, Accuracy) {
    std::vector<double> grid, weight, coef;

    for (int grid_lmax = 17; grid_lmax < 60; grid_lmax +=6) {
        Delley::get(grid_lmax, grid, weight);

        int coef_lmax = grid_lmax / 2;
        randgen(coef_lmax, coef);

        double val = 0.0;
        std::vector<double> ylm_real;
        for (size_t i = 0; i < weight.size(); i++) {
            double tmp = 0.0;
            ModuleBase::Ylm::sph_harm(coef_lmax, grid[3*i], grid[3*i+1], grid[3*i+2], ylm_real);
            for (int l = 0; l <= coef_lmax; l++) {
                for (int m = 0; m <= 2*l; ++m) {
                    tmp += coef[l] * ylm_real[l*l+m];
                }
            }
            val += weight[i] * tmp * tmp;
        }
        val *= 4.0 * std::acos(-1.0);

        double val_ref = 0.0;
        for (int l = 0; l <= coef_lmax; l++) {
            val_ref += (2*l+1) * coef[l] * coef[l];
        }

        double abs_diff = std::abs(val - val_ref);

        printf("order = %2i    abs_diff = %8.5e\n", grid_lmax, abs_diff);
        EXPECT_LT(abs_diff, abs_tol + rel_tol * val_ref);
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
