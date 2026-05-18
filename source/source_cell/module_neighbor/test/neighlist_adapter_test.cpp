#include "../neighlist_adapter.h"

#include "gtest/gtest.h"

#include "prepare_unitcell.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <memory>
#include <vector>

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
    this->start_mag = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_mag;
}

namespace
{

struct AdjEntry
{
    int ntype;
    int natom;
    ModuleBase::Vector3<int> box;
    ModuleBase::Vector3<double> tau;
};

bool AdjEntryLess(const AdjEntry& lhs, const AdjEntry& rhs)
{
    if (lhs.ntype != rhs.ntype)
    {
        return lhs.ntype < rhs.ntype;
    }
    if (lhs.natom != rhs.natom)
    {
        return lhs.natom < rhs.natom;
    }
    if (lhs.box.x != rhs.box.x)
    {
        return lhs.box.x < rhs.box.x;
    }
    if (lhs.box.y != rhs.box.y)
    {
        return lhs.box.y < rhs.box.y;
    }
    return lhs.box.z < rhs.box.z;
}

std::vector<AdjEntry> CollectNonSelfEntries(const AdjacentAtomInfo& adjs)
{
    std::vector<AdjEntry> entries;
    for (int i = 0; i < adjs.adj_num; ++i)
    {
        entries.push_back({adjs.ntype[i], adjs.natom[i], adjs.box[i], adjs.adjacent_tau[i]});
    }
    std::sort(entries.begin(), entries.end(), AdjEntryLess);
    return entries;
}

void ExpectAdjInfoShape(const AdjacentAtomInfo& adjs)
{
    const size_t expected_size = static_cast<size_t>(adjs.adj_num + 1);
    ASSERT_EQ(adjs.ntype.size(), expected_size);
    ASSERT_EQ(adjs.natom.size(), expected_size);
    ASSERT_EQ(adjs.box.size(), expected_size);
    ASSERT_EQ(adjs.adjacent_tau.size(), expected_size);
}

void ExpectSelfLast(const AdjacentAtomInfo& adjs, const UnitCell& ucell, int ntype, int nnumber)
{
    ExpectAdjInfoShape(adjs);
    const int self = adjs.adj_num;
    EXPECT_EQ(adjs.ntype[self], ntype);
    EXPECT_EQ(adjs.natom[self], nnumber);
    EXPECT_EQ(adjs.box[self].x, 0);
    EXPECT_EQ(adjs.box[self].y, 0);
    EXPECT_EQ(adjs.box[self].z, 0);
    EXPECT_NEAR(adjs.adjacent_tau[self].x, ucell.atoms[ntype].tau[nnumber].x, 1e-12);
    EXPECT_NEAR(adjs.adjacent_tau[self].y, ucell.atoms[ntype].tau[nnumber].y, 1e-12);
    EXPECT_NEAR(adjs.adjacent_tau[self].z, ucell.atoms[ntype].tau[nnumber].z, 1e-12);
}

void ExpectEquivalentAdjs(const AdjacentAtomInfo& expected, const AdjacentAtomInfo& actual)
{
    ExpectAdjInfoShape(expected);
    ExpectAdjInfoShape(actual);
    ASSERT_EQ(actual.adj_num, expected.adj_num);

    const std::vector<AdjEntry> expected_entries = CollectNonSelfEntries(expected);
    const std::vector<AdjEntry> actual_entries = CollectNonSelfEntries(actual);
    ASSERT_EQ(actual_entries.size(), expected_entries.size());

    for (size_t i = 0; i < expected_entries.size(); ++i)
    {
        EXPECT_EQ(actual_entries[i].ntype, expected_entries[i].ntype);
        EXPECT_EQ(actual_entries[i].natom, expected_entries[i].natom);
        EXPECT_EQ(actual_entries[i].box.x, expected_entries[i].box.x);
        EXPECT_EQ(actual_entries[i].box.y, expected_entries[i].box.y);
        EXPECT_EQ(actual_entries[i].box.z, expected_entries[i].box.z);
        EXPECT_NEAR(actual_entries[i].tau.x, expected_entries[i].tau.x, 1e-12);
        EXPECT_NEAR(actual_entries[i].tau.y, expected_entries[i].tau.y, 1e-12);
        EXPECT_NEAR(actual_entries[i].tau.z, expected_entries[i].tau.z, 1e-12);
    }
}

AdjacentAtomInfo BuildGridAdjs(const UnitCell& ucell, double radius_lat0, int ntype, int nnumber)
{
    Grid_Driver grid_d(0, 0);
    std::ofstream ofs("neighlist_adapter_test.out");
    grid_d.init(ofs, ucell, radius_lat0, true);
    ofs.close();
    std::remove("neighlist_adapter_test.out");

    AdjacentAtomInfo adjs;
    grid_d.Find_atom(ucell, ntype, nnumber, &adjs);
    return adjs;
}

AdjacentAtomInfo BuildAdapterAdjs(const UnitCell& ucell, double radius_lat0, int ntype, int nnumber)
{
    NeighborSearch neighbor_search;
    neighbor_search.init(ucell, radius_lat0 * ucell.lat0, 0);
    neighbor_search.build_neighbors();

    AdjacentAtomInfo adjs;
    convert_neighbor_search_to_adjs(ucell, neighbor_search, ntype, nnumber, adjs);
    return adjs;
}

void ExpectAdapterMatchesGrid(const UnitCell& ucell, double radius_lat0, int ntype, int nnumber)
{
    const AdjacentAtomInfo expected = BuildGridAdjs(ucell, radius_lat0, ntype, nnumber);
    const AdjacentAtomInfo actual = BuildAdapterAdjs(ucell, radius_lat0, ntype, nnumber);

    ExpectEquivalentAdjs(expected, actual);
    ExpectSelfLast(actual, ucell, ntype, nnumber);
}

bool HasBoxNeighbor(const AdjacentAtomInfo& adjs, int ntype, int nnumber, int cell_x, int cell_y, int cell_z)
{
    for (int i = 0; i < adjs.adj_num; ++i)
    {
        if (adjs.ntype[i] == ntype && adjs.natom[i] == nnumber && adjs.box[i].x == cell_x
            && adjs.box[i].y == cell_y && adjs.box[i].z == cell_z)
        {
            return true;
        }
    }
    return false;
}

} // namespace

TEST(NeighlistAdapterTest, MatchesGridDriverForSiFixture)
{
    std::unique_ptr<UnitCell> ucell(UcellTestLib["Si"].SetUcellInfo());

    ExpectAdapterMatchesGrid(*ucell, 0.5, 0, 0);
    ExpectAdapterMatchesGrid(*ucell, 0.5, 0, 1);
}

TEST(NeighlistAdapterTest, PreservesOrthogonalBoundaryImageBoxes)
{
    UcellTestPrepare prepare("sc",
                             2,
                             false,
                             false,
                             false,
                             "volume",
                             1.0,
                             {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                             {"X"},
                             {"X.upf"},
                             {"upf201"},
                             {"X.orb"},
                             {2},
                             {1.0},
                             "Cartesian",
                             {0.05, 0.0, 0.0, 0.95, 0.0, 0.0});
    std::unique_ptr<UnitCell> ucell(prepare.SetUcellInfo());

    ExpectAdapterMatchesGrid(*ucell, 0.2, 0, 0);
    ExpectAdapterMatchesGrid(*ucell, 0.2, 0, 1);

    const AdjacentAtomInfo atom0_adjs = BuildAdapterAdjs(*ucell, 0.2, 0, 0);
    const AdjacentAtomInfo atom1_adjs = BuildAdapterAdjs(*ucell, 0.2, 0, 1);
    EXPECT_TRUE(HasBoxNeighbor(atom0_adjs, 0, 1, -1, 0, 0));
    EXPECT_TRUE(HasBoxNeighbor(atom1_adjs, 0, 0, 1, 0, 0));
}

TEST(NeighlistAdapterTest, MatchesGridDriverForSkewedCell)
{
    UcellTestPrepare prepare("skew",
                             2,
                             false,
                             false,
                             false,
                             "volume",
                             1.0,
                             {1.0, 0.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
                             {"X"},
                             {"X.upf"},
                             {"upf201"},
                             {"X.orb"},
                             {2},
                             {1.0},
                             "Cartesian",
                             {0.05, 0.0, 0.0, 0.95, 0.0, 0.0});
    std::unique_ptr<UnitCell> ucell(prepare.SetUcellInfo());

    ExpectAdapterMatchesGrid(*ucell, 0.3, 0, 0);
    ExpectAdapterMatchesGrid(*ucell, 0.3, 0, 1);
}
