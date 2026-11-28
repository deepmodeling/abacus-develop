#include "../lri_cv.h"
#include "source_base/parallel_global.h"

#include <RI/physics/Exx.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <utility>
#include <valarray>

namespace
{
template<typename T>
void check_cache_handoff()
{
    using Cell = std::array<int, 3>;
    using Pair = std::pair<int, Cell>;
    using TensorMap = std::map<int, std::map<Pair, RI::Tensor<T>>>;
    LRI_CV<T> cache;
    RI::Exx<int, int, 3, T> engine;
    const std::map<int, std::array<double, 3>> positions{{0, {0.0, 0.0, 0.0}}};
    const std::array<std::array<double, 3>, 3> lattice{{{1.0, 0.0, 0.0},
                                                      {0.0, 1.0, 0.0},
                                                      {0.0, 0.0, 1.0}}};
    const Cell period{{1, 1, 1}};
    engine.set_parallel(MPI_COMM_WORLD, positions, lattice, period);

    // Reuse the same CV object after clearing, as on a later ionic step.
    for (int generation = 1; generation <= 2; ++generation)
    {
        std::weak_ptr<std::valarray<T>> cache_only;
        {
            RI::Tensor<T> c({1, 1, 1});
            RI::Tensor<T> v({1, 1});
            RI::Tensor<T> dc({1, 1, 1});
            RI::Tensor<T> dv({1, 1});
            RI::Tensor<T> unused({2, 2});
            const double coefficient = 0.5 * generation;
            c.ptr()[0] = T(coefficient);
            v.ptr()[0] = T(0.75);
            dc.ptr()[0] = T(0.125);
            dv.ptr()[0] = T(0.25);
            cache.Cws[0][0][{0.0, 0.0, 0.0}] = c;
            cache.Vws[0][0][{0.0, 0.0, 0.0}] = v;
            cache.dCws[0][0][{0.0, 0.0, 0.0}] = {{dc, dc, dc}};
            cache.dVws[0][0][{0.0, 0.0, 0.0}] = {{dv, dv, dv}};
            cache.Vws[1][1][{0.0, 0.0, 0.0}] = unused;
            cache_only = unused.data;

            // Periodic images exercise folding before the retained tensors are used.
            const Pair origin{0, {0, 0, 0}};
            const Pair image{0, {1, 0, 0}};
            const TensorMap cs{{0, {{origin, c}, {image, c}}}};
            const TensorMap vs{{0, {{origin, v}, {image, v}}}};
            const TensorMap dcs{{0, {{origin, dc}, {image, dc}}}};
            const TensorMap dvs{{0, {{origin, dv}, {image, dv}}}};
            const std::array<TensorMap, 3> dc_axes{{dcs, dcs, dcs}};
            const std::array<TensorMap, 3> dv_axes{{dvs, dvs, dvs}};
            engine.set_Cs(cs, 0.0);
            engine.set_Vs(vs, 0.0);
            engine.set_dCs(dc_axes, 0.0);
            engine.set_dVs(dv_axes, 0.0);
            const std::array<std::array<TensorMap, 3>, 3> dcrs{{dc_axes, dc_axes, dc_axes}};
            const std::array<std::array<TensorMap, 3>, 3> dvrs{{dv_axes, dv_axes, dv_axes}};
            engine.set_dCRs(dcrs, 0.0);
            engine.set_dVRs(dvrs, 0.0);
            engine.set_Ds(vs, 0.0);
        }

        engine.cal_Hs();
        engine.cal_force();
        engine.cal_stress();
        const T energy_before = engine.energy;
        const std::array<std::map<int, T>, 3> force_before = engine.force;
        std::array<T, 9> stress_before;
        for (int component = 0; component < 9; ++component)
        {
            stress_before[component] = engine.stress.ptr()[component];
        }
        EXPECT_FALSE(cache_only.expired());
        cache.clear_cache();
        EXPECT_TRUE(cache_only.expired());
        EXPECT_TRUE(cache.Cws.empty());
        EXPECT_TRUE(cache.Vws.empty());
        EXPECT_TRUE(cache.dCws.empty());
        EXPECT_TRUE(cache.dVws.empty());
        cache.clear_cache();

        // Input maps have left scope: these contractions depend on LibRI ownership.
        engine.cal_Hs();
        engine.cal_force();
        engine.cal_stress();
        const T energy_difference = engine.energy - energy_before;
        const double energy_error = std::abs(energy_difference);
        const double energy_tolerance = 1e-12 * (1.0 + std::abs(energy_before));
        EXPECT_LE(energy_error, energy_tolerance);
        EXPECT_EQ(engine.force, force_before);
        for (int component = 0; component < 9; ++component)
        {
            EXPECT_EQ(engine.stress.ptr()[component], stress_before[component]);
        }
        const double energy_magnitude = std::abs(engine.energy);
        EXPECT_GT(energy_magnitude, 0.0);
    }
}
}

TEST(LriCv, ReleaseRealCacheAfterHandoff)
{
    check_cache_handoff<double>();
}

TEST(LriCv, ReleaseComplexCacheAfterHandoff)
{
    check_cache_handoff<std::complex<double>>();
}

int main(int argc, char** argv)
{
    int processes = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, processes, threads, rank);
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
