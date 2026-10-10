#include "../exx_lr_ws.h"
#include "source_base/parallel_global.h"
#include <RI/physics/Exx.h>
#include <gtest/gtest.h>

namespace
{
template<typename T>
void check_workspace()
{
    using Engine = RI::Exx<int, int, 3, T>;
    using Cell = std::array<int, 3>;
    using Pair = std::pair<int, Cell>;
    using TensorMap = std::map<int, std::map<Pair, RI::Tensor<T>>>;
    Engine ks;
    Engine lr;
    const std::map<int, std::array<double, 3>> positions{{0, {0.0, 0.0, 0.0}}};
    const std::array<std::array<double, 3>, 3> lattice{{{1.0, 0.0, 0.0},
                                                     {0.0, 1.0, 0.0},
                                                     {0.0, 0.0, 1.0}}};
    const Cell period{1, 1, 1};
    ks.set_parallel(MPI_COMM_WORLD, positions, lattice, period);
    lr.set_parallel(MPI_COMM_WORLD, positions, lattice, period);
    RI::Tensor<T> c({1, 1, 1});
    RI::Tensor<T> v({1, 1});
    c.ptr()[0] = T(2.0);
    v.ptr()[0] = T(3.0);
    const Pair pair{0, {0, 0, 0}};
    const TensorMap cs{{0, {{pair, c}}}};
    const TensorMap vs{{0, {{pair, v}}}};
    ks.set_Cs(cs, 0.0);
    ks.set_Vs(vs, 0.0);
    const std::array<TensorMap, 3> dcs{{cs, cs, cs}};
    const std::array<TensorMap, 3> dvs{{vs, vs, vs}};
    ks.set_dCs(dcs, 0.0);
    ks.set_dVs(dvs, 0.0);
    const std::array<std::array<TensorMap, 3>, 3> dcrs{{dcs, dcs, dcs}};
    const std::array<std::array<TensorMap, 3>, 3> dvrs{{dvs, dvs, dvs}};
    ks.set_dCRs(dcrs, 0.0);
    ks.set_dVRs(dvrs, 0.0);
    ks.set_Ds(vs, 0.0, "0");
    ks.set_Ds(vs, 0.0, "1");
    const std::array<std::string, 3> gs_names{{"", "", "1"}};
    ks.cal_Hs(gs_names);
    TensorMap gs_hs;
    for (const auto& atom : ks.Hs)
    {
        for (const auto& entry : atom.second) { gs_hs[atom.first][entry.first] = entry.second.copy(); }
    }
    std::map<std::string, TensorMap> original_packs;
    for (const auto& pack : ks.lri.data_pool)
    {
        for (const auto& atom : pack.second.Ds_ab)
        {
            for (const auto& entry : atom.second)
            {
                original_packs[pack.first][atom.first][entry.first] = entry.second.copy();
            }
        }
    }
    const auto gs_bindings = ks.lri.data_ab_name;
    const auto gs_pool_size = ks.lri.data_pool.size();
    const auto gs_energy = ks.energy;
    share_exx_geometry(ks, lr);
    EXPECT_NE(ks.lri.parallel.get(), lr.lri.parallel.get());
    EXPECT_NE(ks.lri.filter_atom.get(), lr.lri.filter_atom.get());
    EXPECT_FALSE(lr.flag_finish.Ds);
    EXPECT_TRUE(lr.Hs.empty());
    EXPECT_EQ(lr.lri.data_pool.count("Ds_0"), 0);
    EXPECT_EQ(lr.lri.data_pool.count("Ds_1"), 0);
    EXPECT_EQ(lr.lri.data_pool.size(), 26);
    for (const auto& pack : lr.lri.data_pool)
    {
        for (const auto& atom : pack.second.Ds_ab)
        {
            for (const auto& entry : atom.second)
            {
                const auto& original = ks.lri.data_pool.at(pack.first).Ds_ab.at(atom.first).at(entry.first);
                EXPECT_EQ(original.ptr(), entry.second.ptr());
            }
        }
    }
    // Exercise the operations that previously replaced the KS density and Hamiltonian.
    RI::Tensor<T> response({1, 1});
    response.ptr()[0] = T(0.5);
    const TensorMap dx{{0, {{pair, response}}}};
    lr.set_Ds(dx, 0.0);
    lr.cal_Hs();
    lr.cal_force();
    EXPECT_EQ(ks.lri.data_pool.size(), gs_pool_size);
    EXPECT_EQ(ks.lri.data_ab_name, gs_bindings);
    EXPECT_EQ(ks.energy, gs_energy);
    for (const auto& pack : original_packs)
    {
        for (const auto& atom : pack.second)
        {
            for (const auto& entry : atom.second)
            {
                const auto& unchanged = ks.lri.data_pool.at(pack.first).Ds_ab.at(atom.first).at(entry.first);
                EXPECT_EQ(unchanged.ptr()[0], entry.second.ptr()[0]);
            }
        }
    }
    EXPECT_EQ(ks.Hs.size(), gs_hs.size());
    for (const auto& atom : gs_hs)
    {
        for (const auto& entry : atom.second)
        {
            const auto& unchanged = ks.Hs.at(atom.first).at(entry.first);
            EXPECT_EQ(unchanged.ptr()[0], entry.second.ptr()[0]);
        }
    }
    EXPECT_EQ(c.ptr()[0], T(2.0));
    EXPECT_EQ(v.ptr()[0], T(3.0));
    // A new ionic-step workspace must see the newly prepared geometry, not old packs.
    RI::Tensor<T> next_c({1, 1, 1});
    next_c.ptr()[0] = T(4.0);
    const TensorMap next_cs{{0, {{pair, next_c}}}};
    ks.set_Cs(next_cs, 0.0);
    Engine next_lr;
    next_lr.set_parallel(MPI_COMM_WORLD, positions, lattice, period);
    share_exx_geometry(ks, next_lr);
    for (const auto& atom : next_lr.lri.data_pool.at("Cs_").Ds_ab)
    {
        for (const auto& entry : atom.second)
        {
            const auto& updated = ks.lri.data_pool.at("Cs_").Ds_ab.at(atom.first).at(entry.first);
            const auto& previous = lr.lri.data_pool.at("Cs_").Ds_ab.at(atom.first).at(entry.first);
            EXPECT_EQ(entry.second.ptr()[0], updated.ptr()[0]);
            EXPECT_NE(entry.second.ptr()[0], previous.ptr()[0]);
        }
    }
}
}

TEST(ExxLRWorkspace, RealGeometryAndElectronicIsolation) { check_workspace<double>(); }
TEST(ExxLRWorkspace, ComplexGeometryAndElectronicIsolation) { check_workspace<std::complex<double>>(); }

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
