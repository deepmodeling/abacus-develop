#include <gtest/gtest.h>

#include "source_cell/md_cell.h"
#include "../neighbor_search.h"

#include <mpi.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace
{
void ensure_mpi_initialized()
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        int provided = 0;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
    }
}

double cell_volume(const ModuleBase::Matrix3& latvec)
{
    const double cx = latvec.e22 * latvec.e33 - latvec.e23 * latvec.e32;
    const double cy = latvec.e23 * latvec.e31 - latvec.e21 * latvec.e33;
    const double cz = latvec.e21 * latvec.e32 - latvec.e22 * latvec.e31;
    return std::abs(latvec.e11 * cx + latvec.e12 * cy + latvec.e13 * cz);
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

MDCell make_mdcell(const ModuleBase::Matrix3& latvec,
                   const std::vector<ModuleBase::Vector3<double> >& positions,
                   double cutoff)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_SELF, &rank);

    const ModuleBase::Matrix3 gt = latvec.Inverse();
    std::vector<LocalAtom> owned_atoms;
    owned_atoms.reserve(positions.size());
    for (std::size_t iat = 0; iat < positions.size(); ++iat)
    {
        owned_atoms.push_back(LocalAtom(positions[iat],
                                        positions[iat] * gt,
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(1, 1, 1),
                                        1.0,
                                        0,
                                        static_cast<int>(iat),
                                        rank,
                                        false));
    }
    return MDCell(latvec,
                  gt,
                  1.0,
                  cell_volume(latvec),
                  static_cast<int>(positions.size()),
                  owned_atoms,
                  std::vector<std::string>(1, "X"),
                  std::vector<double>(1, 1.0),
                  MPI_COMM_SELF,
                  cutoff,
                  0.0);
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

TEST(NeighborSearchTest, MdCellTwoAtomsNeighbor)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};
    MDCell mdcell = make_mdcell(latvec, positions, 1.0);

    NeighborSearch ns;
    ns.init(mdcell, 1.0);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 2);
    EXPECT_EQ(list.get_numneigh(0), 8);
    EXPECT_EQ(list.get_numneigh(1), 8);
}

TEST(NeighborSearchTest, MdCellNoNeighbor)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.0, 0.0, 0.0}, {0.49, 0.0, 0.0}};
    MDCell mdcell = make_mdcell(latvec, positions, 0.1);

    NeighborSearch ns;
    ns.init(mdcell, 0.1);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    ASSERT_EQ(list.get_nlocal(), 2);
    EXPECT_EQ(list.get_numneigh(0), 0);
    EXPECT_EQ(list.get_numneigh(1), 0);
}

TEST(NeighborSearchTest, MdCellInitBuildsOwnedAndGhostAtoms)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};
    MDCell mdcell = make_mdcell(latvec, positions, 1.0);

    NeighborSearch ns;
    ns.init(mdcell, 1.0);

    EXPECT_EQ(mdcell.mpi_size(), 1);
    EXPECT_EQ(ns.get_inside_atoms().size(), 2U);
    EXPECT_GT(ns.get_ghost_atoms().size(), 0U);
    EXPECT_GT(ns.get_all_atoms().size(), ns.get_inside_atoms().size());
    EXPECT_EQ(ns.get_neighbor_list().get_nlocal(), 2);
}

TEST(NeighborSearchTest, MdCellNeighborIdsStayLocalToAllAtoms)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.0, 0.0, 0.0},
                                                              {0.5, 0.0, 0.0},
                                                              {0.0, 0.5, 0.0}};
    MDCell mdcell = make_mdcell(latvec, positions, 0.75);

    NeighborSearch ns;
    ns.init(mdcell, 0.75);
    ns.build_neighbors();

    const NeighborList& list = ns.get_neighbor_list();
    const std::vector<NeighborAtom>& all_atoms = ns.get_all_atoms();
    EXPECT_GT(count_pairs(list), 0U);
    for (std::size_t atom_id = 0; atom_id < all_atoms.size(); ++atom_id)
    {
        EXPECT_EQ(all_atoms[atom_id].atom_id,
                  ModuleNeighList::checked_local_atom_index(atom_id, "test atom id"));
    }
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

TEST(NeighborSearchTest, MdCellPreservesMdAtomStateAcrossOwnedAndGhostStorage)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};
    MDCell mdcell = make_mdcell(latvec, positions, 1.0);

    ASSERT_EQ(mdcell.nlocal(), 2);
    std::vector<LocalAtom>& owned_atoms = mdcell.mutable_owned_atoms();
    owned_atoms[0].vel.set(1.0, 2.0, 3.0);
    owned_atoms[0].mbl.set(1, 0, 1);
    owned_atoms[0].mass = 7.5;
    owned_atoms[1].vel.set(-1.0, -2.0, -3.0);
    owned_atoms[1].mbl.set(0, 1, 1);
    owned_atoms[1].mass = 8.5;

    mdcell.exchange_ghost_atoms();

    ASSERT_GT(mdcell.nghost(), 0);
    const std::vector<LocalAtom>& ghost_atoms = mdcell.ghost_atoms();
    EXPECT_DOUBLE_EQ(ghost_atoms[0].force.x, 0.0);
    EXPECT_DOUBLE_EQ(ghost_atoms[0].force.y, 0.0);
    EXPECT_DOUBLE_EQ(ghost_atoms[0].force.z, 0.0);

    bool found_first = false;
    bool found_second = false;
    for (std::size_t i = 0; i < ghost_atoms.size(); ++i)
    {
        if (ghost_atoms[i].type == owned_atoms[0].type &&
            ghost_atoms[i].type_index == owned_atoms[0].type_index)
        {
            found_first = true;
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.x, 1.0);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.y, 2.0);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.z, 3.0);
            EXPECT_EQ(ghost_atoms[i].mbl.x, 1);
            EXPECT_EQ(ghost_atoms[i].mbl.y, 0);
            EXPECT_EQ(ghost_atoms[i].mbl.z, 1);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].mass, 7.5);
        }
        if (ghost_atoms[i].type == owned_atoms[1].type &&
            ghost_atoms[i].type_index == owned_atoms[1].type_index)
        {
            found_second = true;
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.x, -1.0);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.y, -2.0);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].vel.z, -3.0);
            EXPECT_EQ(ghost_atoms[i].mbl.x, 0);
            EXPECT_EQ(ghost_atoms[i].mbl.y, 1);
            EXPECT_EQ(ghost_atoms[i].mbl.z, 1);
            EXPECT_DOUBLE_EQ(ghost_atoms[i].mass, 8.5);
        }
    }

    EXPECT_TRUE(found_first);
    EXPECT_TRUE(found_second);
}

TEST(NeighborSearchTest, MdCellMigrateOwnedAtomsReassignsOwnership)
{
    ensure_mpi_initialized();

    const ModuleBase::Matrix3 latvec = identity_lattice();
    const std::vector<ModuleBase::Vector3<double> > positions{{0.1, 0.1, 0.1}};
    MDCell mdcell = make_mdcell(latvec, positions, 0.2);

    ASSERT_EQ(mdcell.nlocal(), 1);
    mdcell.mutable_owned_atoms()[0].cart.set(1.2, -0.1, 0.1);
    mdcell.migrate_owned_atoms();

    ASSERT_EQ(mdcell.nlocal(), 1);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].frac.x, 0.2);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].frac.y, 0.9);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].frac.z, 0.1);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].cart.x, 0.2);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].cart.y, 0.9);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].cart.z, 0.1);
}
