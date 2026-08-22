#include <gtest/gtest.h>

#include "source_cell/module_neighlist/local_atom.h"
#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_cell/unitcell.h"

#include <cstddef>
#include <vector>

namespace
{
void initialize_test_ucell(UnitCell& ucell,
                           double lat0,
                           double omega,
                           const ModuleBase::Matrix3& latvec,
                           int ntype,
                           const std::vector<int>& na,
                           const std::vector<ModuleBase::Vector3<double>>& tau)
{
    ucell.lat0 = lat0;
    ucell.omega = omega;
    ucell.latvec = latvec;
    ucell.GT = latvec.Inverse();
    ucell.ntype = ntype;
    ucell.nat = 0;
    ucell.atoms = new Atom[ntype];
    std::size_t iat = 0;
    for (int it = 0; it < ntype; ++it)
    {
        ucell.atoms[it].type = it;
        ucell.atoms[it].na = na[static_cast<std::size_t>(it)];
        ucell.atoms[it].tau.resize(static_cast<std::size_t>(ucell.atoms[it].na));
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ucell.atoms[it].tau[static_cast<std::size_t>(ia)] = tau[iat++];
        }
        ucell.nat += ucell.atoms[it].na;
    }
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

} // namespace

TEST(NeighborSearchTest, TwoAtomsNeighbor)
{
    UnitCell ucell;
    initialize_test_ucell(ucell,
                          1.0,
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
    UnitCell ucell;
    initialize_test_ucell(ucell,
                          1.0,
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
    UnitCell ucell;
    initialize_test_ucell(ucell,
                          1.0,
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
