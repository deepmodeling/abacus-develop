#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>

#include "source_cell/read_stru.h"
#include "source_cell/unitcell.h"

/************************************************
 *  unit test of read_atoms.cpp malformed-input handling
 ***********************************************/

/**
 * - Tested Functions:
 *   - NonNumericCoordinateRejected
 *     - read_atom_positions(): a coordinate line with non-numeric content
 *       (e.g. "1.0 1.0 x") is now rejected at the extraction point; a line
 *       truncated at EOF still goes to the existing mismatch-in-atom-number
 *       exit in parse_atom_properties() (issue #8356)
 *   - NonNumericAtomCountRejected
 *     - read_atom_type_header(): a non-numeric atom count left na = 0 and
 *       triggered the "atom number is 0" banner for the wrong reason; it is
 *       now rejected at the extraction point
 *   - ValidPositionsAccepted
 *     - a well-formed ATOMIC_POSITIONS section still reads normally
 */

class ReadAtomsMalformedTest : public ::testing::Test
{
  protected:
    UnitCell ucell;
    std::ofstream ofs_running;
    std::ofstream ofs_warning;
    std::string log_path;
    std::string warn_path;
    std::string stru_path;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_atoms_malformed.log";
        warn_path = testing::TempDir() + "read_atoms_malformed.warn";
        ofs_running.open(log_path);
        ofs_warning.open(warn_path);
        ucell.ntype = 1;
        ucell.atoms = new Atom[1];
        ucell.atoms[0].label = "H";
        ucell.pseudo_fn.resize(1);
        ucell.pseudo_type.resize(1);
        ucell.orbital_fn.resize(1);
        // a sane lattice scale: finalize_atom_positions() checks the
        // transformed coordinates, which are divided by lat0
        ucell.lat0 = 4.0;
    }

    void TearDown() override
    {
        ofs_running.close();
        ofs_warning.close();
        std::remove(log_path.c_str());
        std::remove(warn_path.c_str());
        if (!stru_path.empty())
        {
            std::remove(stru_path.c_str());
        }
    }

    bool read_positions(const std::string& content)
    {
        stru_path = testing::TempDir() + "positions.stru";
        std::ofstream ofs(stru_path);
        ofs << content;
        ofs.close();
        ofs_running << "created test STRU file: " << stru_path << std::endl;
        std::ifstream ifpos(stru_path);
        const bool result = unitcell::read_atom_positions(ucell,
                                                          ifpos,
                                                          ofs_running,
                                                          ofs_warning,
                                                          1,
                                                          "pw",
                                                          "",
                                                          "atomic",
                                                          0.0,
                                                          false,
                                                          false,
                                                          "scf",
                                                          "ksdft",
                                                          0);
        ifpos.close();
        return result;
    }
};

TEST_F(ReadAtomsMalformedTest, NonNumericCoordinateRejected)
{
    const bool result = read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n2\n0.0 0.0 0.0\n1.0 1.0 x\n");
    EXPECT_FALSE(result);
}

TEST_F(ReadAtomsMalformedTest, NonNumericAtomCountRejected)
{
    const bool result = read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\nx\n0.0 0.0 0.0\n");
    EXPECT_FALSE(result);
}

TEST_F(ReadAtomsMalformedTest, ValidPositionsAccepted)
{
    // the second atom sits at half a cell on each axis: anything on an exact
    // lattice period (e.g. 1.0 1.0 1.0) wraps onto the first atom by PBC and
    // check_tau would legitimately report overlapping atoms
    const bool result =
        read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n2\n0.0 0.0 0.0\n0.5 0.5 0.5\n");
    EXPECT_TRUE(result);
    EXPECT_EQ(ucell.atoms[0].na, 2);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].tau[1].z, 0.5);
}
