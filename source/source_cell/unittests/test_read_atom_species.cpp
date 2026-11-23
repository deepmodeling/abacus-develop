#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "source_cell/read_stru.h"
#include "source_cell/unitcell.h"

/************************************************
 *  unit test of read_atom_species.cpp malformed-input handling
 ***********************************************/

/**
 * - Tested Functions:
 *   - SpeciesMissingMassRejected
 *     - read_atom_species(): an ATOMIC_SPECIES line without a mass left
 *       mass = 0 silently (issue #8356); it must now quit with a message
 *   - PartialLatticeVectorsRejected
 *     - read_lattice_constant(): fewer than 9 LATTICE_VECTORS components left
 *       the unread tail at the Matrix3 default (identity), surviving even the
 *       determinant check; it must now quit with a message
 *   - ValidSpeciesAndVectors
 *     - a well-formed species line and a full 3x3 LATTICE_VECTORS block still
 *       read normally
 *   - LatticeVectorsWithInlineCommentsAccepted
 *     - vector rows annotated with trailing comments ("0.5 0.5 0.0  # a3"),
 *       common in real STRU files, must still be accepted
 */

class ReadSpeciesMalformedTest : public ::testing::Test
{
  protected:
    UnitCell ucell;
    std::ofstream ofs_running;
    std::string log_path;
    std::vector<std::string> stru_paths;
    std::unique_ptr<Atom[]> owned_atoms;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_species_malformed.log";
        ofs_running.open(log_path);
        ucell.ntype = 1;
        // UnitCell::~UnitCell() deletes atoms only when set_atom_flag is
        // true; keep ownership here instead of relying on that flag
        owned_atoms.reset(new Atom[1]);
        ucell.atoms = owned_atoms.get();
        ucell.pseudo_fn.resize(1);
        ucell.pseudo_type.resize(1);
        ucell.orbital_fn.resize(1);
    }

    void TearDown() override
    {
        ofs_running.close();
        std::remove(log_path.c_str());
        for (const std::string& path : stru_paths)
        {
            std::remove(path.c_str());
        }
    }

    std::string write_stru_file(const std::string& name, const std::string& content)
    {
        const std::string path = testing::TempDir() + name;
        std::ofstream ofs(path);
        ofs << content;
        ofs.close();
        stru_paths.push_back(path);
        ofs_running << "created test STRU file: " << path << std::endl;
        return path;
    }

    std::ifstream open_stru_file(const std::string& path)
    {
        const std::string msg = "reading test STRU file: " + path;
        std::cout << msg << std::endl;
        ofs_running << msg << std::endl;
        return std::ifstream(path);
    }
};

TEST_F(ReadSpeciesMalformedTest, SpeciesMissingMassRejected)
{
    const std::string stru = write_stru_file("species_no_mass.stru", "ATOMIC_SPECIES\nH\n");
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_atom_species(ifa, ofs_running, ucell, "pw", "", "atomic", 0.0, false, false),
                ::testing::ExitedWithCode(1),
                "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid ATOMIC_SPECIES line"));
    ifa.close();
}

TEST_F(ReadSpeciesMalformedTest, PartialLatticeVectorsRejected)
{
    const std::string stru = write_stru_file("latvec_partial.stru",
                                             "LATTICE_CONSTANT\n4.0\nLATTICE_VECTORS\n1 0 0\n0 1 0\n");
    ucell.latName = "user_defined_lattice";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, ucell.lat),
                ::testing::ExitedWithCode(1),
                "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_VECTORS"));
    ifa.close();
}

TEST_F(ReadSpeciesMalformedTest, ValidSpeciesAndVectors)
{
    const std::string stru = write_stru_file("species_ok.stru", "ATOMIC_SPECIES\nH 1.008 H.upf\n");
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_atom_species(ifa, ofs_running, ucell, "pw", "", "atomic", 0.0, false, false));
    EXPECT_DOUBLE_EQ(ucell.atoms[0].mass, 1.008);
    ifa.close();

    const std::string stru2 = write_stru_file("latvec_ok.stru",
                                              "LATTICE_CONSTANT\n4.0\nLATTICE_VECTORS\n1 0 0\n0 1 0\n0 0 1\n");
    ucell.latName = "user_defined_lattice";
    std::ifstream ifa2 = open_stru_file(stru2);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa2, ofs_running, ucell.lat));
    EXPECT_DOUBLE_EQ(ucell.lat.latvec.e33, 1.0);
    ifa2.close();
}

TEST_F(ReadSpeciesMalformedTest, LatticeVectorsWithInlineCommentsAccepted)
{
    // real STRU files annotate vector rows with trailing comments
    // (e.g. "0.5 0.5 0.0  # a3"); the reader must swallow the rest of
    // each row like the old READ_VALUE-based code did
    const std::string stru = write_stru_file("latvec_comments.stru",
                                             "LATTICE_CONSTANT\n4.0\nLATTICE_VECTORS\n"
                                             "1 0 0  # Lattice vector 1\n"
                                             "0 1 0  # Lattice vector 2\n"
                                             "0 0 1  # Lattice vector 3\n");
    ucell.latName = "user_defined_lattice";
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, ucell.lat));
    EXPECT_DOUBLE_EQ(ucell.lat.latvec.e33, 1.0);
    ifa.close();
}
