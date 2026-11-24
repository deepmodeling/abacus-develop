#include "source_cell/mdcell_adapter.h"

#include "source_base/matrix.h"
#include "source_base/parallel_cell.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_cell/unitcell.h"

#include "gtest/gtest.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

TEST(MDCellAdapterTest, StruRestartDoesNotReplicateCheckpointAgain)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const long long stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const int rank = domain.rank();
    const std::string file_prefix = testing::TempDir() + "mdcell_adapter_" + std::to_string(rank)
                                    + "_" + std::to_string(stamp);
    const std::string stru_file = file_prefix + ".STRU";
    const std::string log_file = file_prefix + ".log";
    const char* const stru_path = stru_file.c_str();
    const char* const log_path = log_file.c_str();
    {
        std::ofstream out(stru_path);
        ASSERT_TRUE(out.good());
        out << "ATOMIC_SPECIES\nHe 4.0026 auto auto\n\n"
               "LATTICE_CONSTANT\n1.0\n\n"
               "LATTICE_VECTORS\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\n\n"
               "ATOMIC_POSITIONS\nDirect\n\nHe\n0.0\n1\n"
               "0.25 0.25 0.25 m 1 1 1\n";
        ASSERT_TRUE(out.good());
    }

    std::ofstream running_log(log_path);
    ASSERT_TRUE(running_log.good());
    const std::vector<int> replica = {2, 1, 1};
    MDCell fresh;
    DomainDecomposition fresh_decomp;
    MDCellAdapter::initialize_from_stru(fresh, fresh_decomp, stru_file, replica, false, 1.0, running_log);
    EXPECT_EQ(fresh.nat(), 2);

    MDCell restarted;
    DomainDecomposition restart_decomp;
    MDCellAdapter::initialize_from_stru(restarted, restart_decomp, stru_file, replica, true, 1.0, running_log);
    EXPECT_EQ(restarted.nat(), 1);
    running_log.close();

    std::ifstream log_input(log_path);
    ASSERT_TRUE(log_input.good());
    std::string log_line;
    int atom_count_lines = 0;
    while (std::getline(log_input, log_line))
    {
        const bool is_atom_count = log_line.find("TOTAL ATOM NUMBER") != std::string::npos;
        if (is_atom_count)
        {
            ++atom_count_lines;
        }
    }
    EXPECT_EQ(atom_count_lines, 2);

    std::remove(stru_path);
    std::remove(log_path);
}

TEST(MDCellAdapterTest, InitializeFromUnitCellPreservesAtomsAndAllowsBacking)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e12 = 0.0;
    ucell.latvec.e13 = 0.0;
    ucell.latvec.e21 = 0.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e23 = 0.0;
    ucell.latvec.e31 = 0.0;
    ucell.latvec.e32 = 0.0;
    ucell.latvec.e33 = 1.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.atoms = new Atom[ucell.ntype];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Ar";
    ucell.atoms[0].mass = 39.948;
    ucell.atoms[0].na = 1;
    ucell.atoms[0].tau.resize(1);
    ucell.atoms[0].taud.resize(1);
    ucell.atoms[0].vel.resize(1);
    ucell.atoms[0].mbl.resize(1);
    ucell.atoms[0].tau[0].set(-0.2, 1.25, 2.5);
    ucell.atoms[0].taud[0] = ucell.atoms[0].tau[0];
    ucell.atoms[0].vel[0].set(1.0, 2.0, 3.0);
    ucell.atoms[0].mbl[0].set(0, 0, 0);

    MDCell mdcell;
    DomainDecomposition decomp;
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    MDCellAdapter::initialize_from_unitcell(mdcell, decomp, ucell, 0.0, domain);

    EXPECT_EQ(mdcell.nat(), ucell.nat);
    EXPECT_EQ(mdcell.stru_meta().species.size(), 1U);
    EXPECT_FALSE(mdcell.has_backing_unitcell());
    MDCellAdapter::initialize_backed_from_unitcell(mdcell, decomp, ucell);
    EXPECT_EQ(&mdcell.backing_unitcell(), &ucell);
    const ModuleBase::Vector3<double> frac(0.8, 0.25, 0.5);
    const int owner = decomp.owner_rank_from_frac(frac);
    if (domain.rank() != owner)
    {
        EXPECT_TRUE(mdcell.owned_atoms().empty());
        return;
    }
    ASSERT_EQ(mdcell.owned_atoms().size(), 1U);
    const LocalAtom& atom = mdcell.owned_atoms()[0];
    EXPECT_NEAR(atom.frac.x, 0.8, 1.0e-12);
    EXPECT_NEAR(atom.frac.y, 0.25, 1.0e-12);
    EXPECT_NEAR(atom.frac.z, 0.5, 1.0e-12);
    EXPECT_DOUBLE_EQ((atom.cart - atom.frac * mdcell.latvec()).norm2(), 0.0);
    EXPECT_DOUBLE_EQ(atom.vel.y, 2.0);
    EXPECT_EQ(atom.mbl.x, 0);
    EXPECT_EQ(atom.owner_rank, owner);
}

TEST(MDCellAdapterTest, SyncGeometryAndCollectForces)
{
    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = 2;
    ucell.lat0 = 1.0;
    ucell.latvec.e11 = 1.0;
    ucell.latvec.e22 = 1.0;
    ucell.latvec.e33 = 1.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.omega = 1.0;
    ucell.atoms = new Atom[1];
    ucell.set_atom_flag = true;
    ucell.atoms[0].label = "Ar";
    ucell.atoms[0].mass = 39.948;
    ucell.atoms[0].na = 2;
    ucell.atoms[0].tau.resize(2);
    ucell.atoms[0].taud.resize(2);
    ucell.atoms[0].vel.resize(2);
    ucell.atoms[0].mbl.resize(2);
    ucell.atoms[0].tau[0].set(0.1, 0.1, 0.1);
    ucell.atoms[0].tau[1].set(0.7, 0.7, 0.7);

    MDCell mdcell;
    DomainDecomposition decomp;
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    MDCellAdapter::initialize_from_unitcell(mdcell, decomp, ucell, 0.0, domain);
    mdcell.set_neighbor_cutoff(0.4);

    ucell.lat0 = 1.5;
    ucell.latvec.e11 = 2.0;
    ucell.latvec.e22 = 2.0;
    ucell.latvec.e33 = 2.0;
    ucell.GT = ucell.latvec.Inverse();
    ucell.omega = 27.0;
    ucell.atoms[0].tau[0].set(0.3, 0.4, 0.5);
    ucell.atoms[0].tau[1].set(1.4, 1.5, 1.6);
    MDCellAdapter::sync_geometry_from_unitcell(mdcell, ucell);
    EXPECT_DOUBLE_EQ(mdcell.lat0(), 1.5);
    EXPECT_DOUBLE_EQ(mdcell.omega(), 27.0);
    for (const LocalAtom& atom : mdcell.owned_atoms())
    {
        const double expected_x = atom.type_index == 0 ? 0.3 : 1.4;
        EXPECT_NEAR(atom.cart.x, expected_x, 1.0e-12);
    }

    decomp.prepare_neighbors(mdcell);
    for (LocalAtom& atom : mdcell.owned_atoms())
    {
        const double force = static_cast<double>(atom.type_index + 1);
        atom.force.set(force, 2.0 * force, 3.0 * force);
    }
    ModuleBase::matrix full_force;
    MDCellAdapter::collect_full_forces(mdcell, decomp, ucell, full_force);
    EXPECT_DOUBLE_EQ(full_force(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(full_force(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(full_force(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(full_force(1, 2), 6.0);
}
