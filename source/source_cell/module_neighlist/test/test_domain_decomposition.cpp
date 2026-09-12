#include "source_cell/module_neighlist/domain_decomposition.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_cell/mdcell.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_reduce.h"

#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

namespace
{
class DomainDecompositionTest : public testing::Test
{
protected:
    MDCell cell;
    DomainDecomposition decomp;
    ModuleBase::CommunicationDomain domain;

    void SetUp() override
    {
        domain = ModuleBase::world_comm_domain();
        ModuleBase::Matrix3 lattice;
        lattice.e11 = 4.0;
        lattice.e22 = 4.0;
        lattice.e33 = 4.0;
        std::vector<LocalAtom> atoms;
        // Start with all atoms on rank zero; prepare_neighbors must distribute them.
        if (domain.rank() == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                LocalAtom atom;
                atom.frac.set(i == 0 ? 0.05 : 0.95, 0.5, 0.5);
                atom.cart = atom.frac * lattice;
                atom.type_index = i;
                atom.vel.set(i + 1.0, 2.0, 3.0);
                atom.mbl.set(1, 0, 1);
                atoms.push_back(atom);
            }
        }
        cell.initialize_from_owned_atoms(lattice, lattice.Inverse(), 1.0, 64.0, 2, atoms,
                                          std::vector<std::string>(1, "X"),
                                          std::vector<double>(1, 1.0),
                                          std::vector<std::int64_t>(1, 2), 0.4, domain);
        decomp.init(domain, lattice, 1.0, 0.0, 0.4);
        cell.set_neighbor_cutoff(0.5);
    }

    void expect_neighbor_count(int count)
    {
        ASSERT_TRUE(cell.has_neighbor_search());
        const NeighborList& list = cell.neighbor_search().get_neighbor_list();
        ASSERT_EQ(list.get_ncentral_atoms(), cell.owned_atoms().size());
        for (int i = 0; i < cell.owned_atoms().size(); ++i)
        {
            EXPECT_EQ(list.get_numneigh(i), count);
        }
        long long owned = cell.owned_atoms().size();
        Parallel_Reduce::reduce_all(owned);
        EXPECT_EQ(owned, 2);
    }
};

TEST_F(DomainDecompositionTest, CutoffSetterDoesNotExchangeAtoms)
{
    EXPECT_FALSE(cell.has_neighbor_search());
    EXPECT_EQ(cell.ghost_atoms().size(), 0);
    EXPECT_EQ(cell.owned_atoms().size(), domain.rank() == 0 ? 2 : 0);
    EXPECT_THROW(cell.set_neighbor_cutoff(0.0), std::runtime_error);
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(1);
    for (const LocalAtom& atom : cell.owned_atoms())
    {
        EXPECT_EQ(atom.owner_rank, domain.rank());
        EXPECT_DOUBLE_EQ(atom.vel.x, atom.type_index + 1.0);
        EXPECT_EQ(atom.mbl.y, 0);
    }
}

TEST_F(DomainDecompositionTest, GhostForcesAreReturnedAndConsumed)
{
    decomp.prepare_neighbors(cell);
    long long copies[2] = {0, 0};
    for (LocalAtom& atom : cell.ghost_atoms())
    {
        ++copies[atom.type_index];
        atom.force.set(atom.type_index + 1.0, 0.0, 0.0);
    }
    Parallel_Reduce::reduce_all(copies, 2);
    for (LocalAtom& atom : cell.owned_atoms())
    {
        atom.force.set(5.0, 0.0, 0.0);
    }
    decomp.accumulate_ghost_forces(cell);
    for (const LocalAtom& atom : cell.owned_atoms())
    {
        EXPECT_DOUBLE_EQ(atom.force.x, 5.0 + copies[atom.type_index] * (atom.type_index + 1.0));
    }
    for (const LocalAtom& atom : cell.ghost_atoms())
    {
        EXPECT_DOUBLE_EQ(atom.force.norm2(), 0.0);
    }
    decomp.accumulate_ghost_forces(cell);
    for (const LocalAtom& atom : cell.owned_atoms())
    {
        EXPECT_DOUBLE_EQ(atom.force.x, 5.0 + copies[atom.type_index] * (atom.type_index + 1.0));
    }
}

TEST_F(DomainDecompositionTest, ReusesLayoutAndRefreshesGhostCoordinates)
{
    decomp.prepare_neighbors(cell);
    const NeighborSearch* search = &cell.neighbor_search();
    const std::size_t ghosts = cell.ghost_atoms().size();
    for (LocalAtom& atom : cell.owned_atoms())
    {
        atom.frac.x += 0.01;
        atom.cart = atom.frac * cell.latvec();
    }
    for (LocalAtom& atom : cell.ghost_atoms())
    {
        atom.force.set(8.0, 0.0, 0.0);
    }
    decomp.prepare_neighbors(cell);
    EXPECT_EQ(&cell.neighbor_search(), search);
    EXPECT_EQ(cell.ghost_atoms().size(), ghosts);
    for (const LocalAtom& atom : cell.ghost_atoms())
    {
        EXPECT_NEAR(atom.frac.x, atom.type_index == 0 ? 0.06 : 0.96, 1.0e-12);
        EXPECT_DOUBLE_EQ(atom.force.norm2(), 0.0);
    }
    expect_neighbor_count(1);
}

TEST_F(DomainDecompositionTest, LatticeAndCutoffChangesRebuildNeighbors)
{
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(1);
    ModuleBase::Matrix3 lattice = cell.latvec();
    lattice.e11 = 8.0;
    cell.set_lattice_vectors(lattice);
    cell.refresh_cart_from_frac();
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(0);
    cell.set_neighbor_cutoff(1.0);
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(1);
}

TEST_F(DomainDecompositionTest, OneRankInvalidationRebuildsCollectively)
{
    decomp.prepare_neighbors(cell);
    if (domain.rank() == 0)
    {
        cell.refresh_cart_from_frac();
    }
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(1);
}

TEST_F(DomainDecompositionTest, CrossingSkinThresholdChangesActiveNeighbors)
{
    decomp.prepare_neighbors(cell);
    for (LocalAtom& atom : cell.owned_atoms())
    {
        if (atom.type_index == 0)
        {
            atom.frac.x = 0.25;
            atom.cart = atom.frac * cell.latvec();
        }
    }
    decomp.prepare_neighbors(cell);
    expect_neighbor_count(0);
}

TEST_F(DomainDecompositionTest, SingleDomainUsesHaloFiltering)
{
    EXPECT_EQ(decomp.dims()[0] * decomp.dims()[1] * decomp.dims()[2], domain.size());
    if (domain.size() != 1)
    {
        return;
    }
    const std::array<int, 3> dims = {{1, 1, 1}};
    const std::array<int, 3> coords = {{0, 0, 0}};
    EXPECT_EQ(decomp.dims(), dims);
    EXPECT_EQ(decomp.coords(), coords);
    decomp.prepare_neighbors(cell);
    // Only the two x-boundary images overlap the halo, not all 26 images per atom.
    ASSERT_EQ(cell.ghost_atoms().size(), 2);
    for (const LocalAtom& ghost : cell.ghost_atoms())
    {
        EXPECT_EQ(ghost.owner_rank, 0);
        EXPECT_NEAR(ghost.cart.x, ghost.type_index == 0 ? 4.2 : -0.2, 1.0e-12);
    }
    expect_neighbor_count(1);
}

TEST_F(DomainDecompositionTest, CrossingHaloEdgePreservesCachedGhostSlots)
{
    for (LocalAtom& atom : cell.owned_atoms())
    {
        if (atom.type_index == 0)
        {
            atom.frac.x = 0.224; // Just inside the 0.225 fractional halo margin.
            atom.cart = atom.frac * cell.latvec();
        }
    }
    decomp.prepare_neighbors(cell);
    const std::vector<LocalAtom> ghosts = cell.ghost_atoms();
    const NeighborSearch* search = &cell.neighbor_search();
    for (LocalAtom& atom : cell.owned_atoms())
    {
        if (atom.type_index == 0)
        {
            atom.frac.x += 0.002; // Cross the halo edge without consuming half the skin.
            atom.cart = atom.frac * cell.latvec();
        }
    }
    decomp.prepare_neighbors(cell);
    EXPECT_EQ(&cell.neighbor_search(), search);
    ASSERT_EQ(cell.ghost_atoms().size(), ghosts.size());
    for (std::size_t i = 0; i < ghosts.size(); ++i)
    {
        const LocalAtom& updated = cell.ghost_atoms()[i];
        EXPECT_EQ(updated.type_index, ghosts[i].type_index);
        EXPECT_EQ(updated.owner_rank, ghosts[i].owner_rank);
        EXPECT_NEAR(updated.cart.x - ghosts[i].cart.x,
                    updated.type_index == 0 ? 0.008 : 0.0, 1.0e-12);
        EXPECT_NEAR(updated.cart.y, ghosts[i].cart.y, 1.0e-12);
        EXPECT_NEAR(updated.cart.z, ghosts[i].cart.z, 1.0e-12);
    }
}

TEST_F(DomainDecompositionTest, SkewCellMultipleImagesMatchBruteForce)
{
    ModuleBase::Matrix3 lattice = cell.latvec();
    lattice.e21 = 0.8;
    lattice.e31 = 0.4;
    lattice.e32 = 0.6;
    cell.set_lattice_vectors(lattice);
    cell.refresh_cart_from_frac();
    const double cutoff = 4.5; // Larger than a cell height: multiple image layers.
    cell.set_neighbor_cutoff(cutoff);
    for (int evaluation = 0; evaluation < 2; ++evaluation)
    {
        if (evaluation == 1)
        {
            for (LocalAtom& atom : cell.owned_atoms())
            {
                atom.frac.x += 0.01;
                atom.cart = atom.frac * lattice;
            }
        }
        decomp.prepare_neighbors(cell);
        const NeighborList& list = cell.neighbor_search().get_neighbor_list();
        for (int i = 0; i < cell.owned_atoms().size(); ++i)
        {
            const LocalAtom& center = cell.owned_atoms()[i];
            // Identify each neighbor by atom ID and its integer periodic shift.
            std::set<std::array<int, 4>> expected;
            for (int id = 0; id < 2; ++id)
            {
                const ModuleBase::Vector3<double> frac((id == 0 ? 0.05 : 0.95) + evaluation * 0.01,
                                                       0.5, 0.5);
                for (int x = -3; x <= 3; ++x)
                {
                    for (int y = -3; y <= 3; ++y)
                    {
                        for (int z = -3; z <= 3; ++z)
                        {
                            if (id == center.type_index && x == 0 && y == 0 && z == 0)
                            {
                                continue;
                            }
                            const ModuleBase::Vector3<double> image = frac + ModuleBase::Vector3<double>(x, y, z);
                            if (((image - center.frac) * lattice).norm2() < cutoff * cutoff)
                            {
                                expected.insert({{id, x, y, z}});
                            }
                        }
                    }
                }
            }
            std::set<std::array<int, 4>> actual;
            const int* neighbors = list.get_firstneigh(i);
            for (int j = 0; j < list.get_numneigh(i); ++j)
            {
                const int index = neighbors[j];
                const LocalAtom& atom = index < cell.owned_atoms().size()
                                            ? cell.owned_atoms()[index]
                                            : cell.ghost_atoms()[index - cell.owned_atoms().size()];
                const ModuleBase::Vector3<double> shift = atom.cart * lattice.Inverse() - atom.frac;
                actual.insert({{static_cast<int>(atom.type_index),
                                static_cast<int>(std::lround(shift.x)),
                                static_cast<int>(std::lround(shift.y)),
                                static_cast<int>(std::lround(shift.z))}});
            }
            EXPECT_EQ(actual.size(), static_cast<std::size_t>(list.get_numneigh(i)));
            EXPECT_EQ(actual, expected);
        }
    }
}
} // namespace
