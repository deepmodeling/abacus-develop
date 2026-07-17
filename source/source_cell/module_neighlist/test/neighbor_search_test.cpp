#include <gtest/gtest.h>

#include "../local_atom.h"
#include "../neighbor_search.h"
#include "../unitcell_lite.h"

#include <cstddef>
#include <vector>

namespace
{
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
