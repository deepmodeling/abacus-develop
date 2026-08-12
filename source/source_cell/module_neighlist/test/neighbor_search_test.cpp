#include <gtest/gtest.h>

#include "source_cell/module_neighlist/local_atom.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_cell/module_neighlist/unitcell_lite.h"

#include <cstddef>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
class CountingAtomProvider : public AtomProvider
{
public:
    explicit CountingAtomProvider(const UnitCellLite& ucell) : ucell_(ucell) {}

    double get_lat0() const override
    {
        return ucell_.get_lat0();
    }

    double get_omega() const override
    {
        return ucell_.get_omega();
    }

    const ModuleBase::Matrix3& get_latvec() const override
    {
        return ucell_.get_latvec();
    }

    int get_natom() const override
    {
        return ucell_.get_natom();
    }

    int get_na(int type) const override
    {
        return ucell_.get_na(type);
    }

    int get_ntype() const override
    {
        return ucell_.get_ntype();
    }

    ModuleBase::Vector3<double> get_tau(int type, int index) const override
    {
        ++get_tau_calls_;
        return ucell_.get_tau(type, index);
    }

    std::size_t get_tau_calls() const
    {
        return get_tau_calls_;
    }

private:
    const UnitCellLite& ucell_;
    mutable std::size_t get_tau_calls_ = 0;
};

UnitCellLite make_test_ucell(double lat0,
                             double omega,
                             const ModuleBase::Matrix3& latvec,
                             int ntype,
                             const std::vector<int>& na,
                             const std::vector<ModuleBase::Vector3<double>>& tau)
{
    UnitCellLite ucell;
    ucell.set_lattice(lat0, omega, latvec);
    ucell.set_atoms(ntype, na, tau);
    return ucell;
}

ModuleBase::Matrix3 identity_lattice()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1.0;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = 0.0;
    latvec.e22 = 1.0;
    latvec.e23 = 0.0;
    latvec.e31 = 0.0;
    latvec.e32 = 0.0;
    latvec.e33 = 1.0;
    return latvec;
}

std::size_t count_pairs(const NeighborList& list)
{
    std::size_t pairs = 0;
    for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
    {
        pairs += static_cast<std::size_t>(list.get_numneigh(local_i));
    }
    return pairs;
}

#ifdef _OPENMP
void expect_same_atoms(const std::vector<NeighborAtom>& lhs,
                       const std::vector<NeighborAtom>& rhs)
{
    ASSERT_EQ(lhs.size(), rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(lhs[i].position_x, rhs[i].position_x) << "atom " << i;
        EXPECT_DOUBLE_EQ(lhs[i].position_y, rhs[i].position_y) << "atom " << i;
        EXPECT_DOUBLE_EQ(lhs[i].position_z, rhs[i].position_z) << "atom " << i;
        EXPECT_EQ(lhs[i].atom_type, rhs[i].atom_type) << "atom " << i;
        EXPECT_EQ(lhs[i].atom_index, rhs[i].atom_index) << "atom " << i;
        EXPECT_EQ(lhs[i].atom_id, rhs[i].atom_id) << "atom " << i;
        EXPECT_EQ(lhs[i].global_id, rhs[i].global_id) << "atom " << i;
        EXPECT_EQ(lhs[i].owner_rank, rhs[i].owner_rank) << "atom " << i;
    }
}
#endif
} // namespace

TEST(NeighborSearchTest, TwoAtomsNeighbor)
{
    UnitCellLite ucell = make_test_ucell(1.0,
                                         1.0,
                                         identity_lattice(),
                                         1,
                                         {2},
                                         {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}});

    NeighborSearch ns;
    ns.init(ucell, 1.0);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 2);
    EXPECT_EQ(list.get_numneigh(0), 8);
    EXPECT_EQ(list.get_numneigh(1), 8);
}

TEST(NeighborSearchTest, NoNeighbor)
{
    UnitCellLite ucell = make_test_ucell(1.0,
                                         1.0,
                                         identity_lattice(),
                                         1,
                                         {2},
                                         {{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}});

    NeighborSearch ns;
    ns.init(ucell, 0.1);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 2);
    EXPECT_EQ(list.get_numneigh(0), 0);
    EXPECT_EQ(list.get_numneigh(1), 0);
}

TEST(NeighborSearchTest, SerialInitOwnsCentralAtomsAndBuildsImages)
{
    UnitCellLite ucell = make_test_ucell(1.0,
                                         1.0,
                                         identity_lattice(),
                                         1,
                                         {2},
                                         {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}});

    NeighborSearch ns;
    ns.init(ucell, 1.0);

    EXPECT_EQ(ns.get_inside_atoms().size(), 2U);
    EXPECT_EQ(ns.get_neighbor_list().get_nlocal(), 2);
    EXPECT_EQ(ns.get_all_atoms().size(), 54U);

    const std::vector<NeighborAtom>& all_atoms = ns.get_all_atoms();
    for (std::size_t i = 0; i < all_atoms.size(); ++i)
    {
        EXPECT_EQ(all_atoms[i].atom_id,
                  ModuleNeighList::checked_local_atom_index(i, "test atom id"));
    }
}

TEST(NeighborSearchTest, SerialInitReadsEachPrimaryCoordinateOnce)
{
    UnitCellLite ucell = make_test_ucell(1.0,
                                         1.0,
                                         identity_lattice(),
                                         1,
                                         {2},
                                         {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}});
    CountingAtomProvider provider(ucell);

    NeighborSearch ns;
    ns.init(provider, 1.0);

    EXPECT_EQ(provider.get_tau_calls(), 2U);
    EXPECT_EQ(ns.get_inside_atoms().size(), 2U);
    EXPECT_EQ(ns.get_ghost_atoms().size(), 52U);
    EXPECT_EQ(ns.get_all_atoms().size(), 54U);
}

#ifdef _OPENMP
TEST(NeighborSearchTest, ParallelImageGenerationMatchesSerialOrder)
{
    const int atom_count = 4000;
    std::vector<ModuleBase::Vector3<double>> positions;
    positions.reserve(atom_count);
    for (int atom = 0; atom < atom_count; ++atom)
    {
        positions.push_back(ModuleBase::Vector3<double>(0.0001 * atom,
                                                        0.0002 * atom,
                                                        0.0003 * atom));
    }
    UnitCellLite ucell = make_test_ucell(1.0,
                                         1.0,
                                         identity_lattice(),
                                         1,
                                         {atom_count},
                                         positions);

    const int previous_dynamic = omp_get_dynamic();
    const int previous_max_threads = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(1);
    NeighborSearch serial;
    serial.init(ucell, 1.0);

    omp_set_num_threads(4);
    NeighborSearch parallel;
    parallel.init(ucell, 1.0);

    omp_set_num_threads(previous_max_threads);
    omp_set_dynamic(previous_dynamic);

    EXPECT_EQ(parallel.get_inside_atoms().size(), static_cast<std::size_t>(atom_count));
    EXPECT_EQ(parallel.get_ghost_atoms().size(), 104000U);
    expect_same_atoms(serial.get_inside_atoms(), parallel.get_inside_atoms());
    expect_same_atoms(serial.get_ghost_atoms(), parallel.get_ghost_atoms());
    expect_same_atoms(serial.get_all_atoms(), parallel.get_all_atoms());
}
#endif

TEST(NeighborSearchTest, DistributedInputUsesOwnedCentersAndGhostNeighbors)
{
    std::vector<LocalAtom> owned_atoms;
    std::vector<LocalAtom> ghost_atoms;
    owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    0,
                                    0,
                                    0,
                                    0,
                                    false));
    ghost_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    0,
                                    1,
                                    1,
                                    1,
                                    true));

    NeighborSearch ns;
    ns.init_distributed(owned_atoms, ghost_atoms, 1.0, 1.0);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 1);
    ASSERT_EQ(list.get_numneigh(0), 1);

    const int neighbor_id = list.get_firstneigh(0)[0];
    ASSERT_GE(neighbor_id, 0);
    ASSERT_LT(static_cast<std::size_t>(neighbor_id), ns.get_all_atoms().size());
    EXPECT_EQ(ns.get_all_atoms()[neighbor_id].global_id, 1);
    EXPECT_EQ(ns.get_all_atoms()[neighbor_id].owner_rank, 1);
}

TEST(NeighborSearchTest, DistributedNeighborIdsStayLocalToAllAtoms)
{
    std::vector<LocalAtom> owned_atoms;
    std::vector<LocalAtom> ghost_atoms;
    owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    0,
                                    10,
                                    0,
                                    0,
                                    false));
    owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(2.0, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(2.0, 0.0, 0.0),
                                    0,
                                    11,
                                    1,
                                    0,
                                    false));
    ghost_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.5, 0.0, 0.0),
                                    0,
                                    20,
                                    2,
                                    1,
                                    true));

    NeighborSearch ns;
    ns.init_distributed(owned_atoms, ghost_atoms, 0.75, 1.0);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    const std::vector<NeighborAtom>& all_atoms = ns.get_all_atoms();
    EXPECT_EQ(count_pairs(list), 1U);
    for (int local_i = 0; local_i < list.get_nlocal(); ++local_i)
    {
        for (int ad = 0; ad < list.get_numneigh(local_i); ++ad)
        {
            const int neighbor_id = list.get_firstneigh(local_i)[ad];
            EXPECT_GE(neighbor_id, 0);
            EXPECT_LT(static_cast<std::size_t>(neighbor_id), all_atoms.size());
        }
    }
}
