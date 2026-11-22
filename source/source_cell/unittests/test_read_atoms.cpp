#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

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
 *     - read_atom_type_header(): a non-numeric atom count after a valid first
 *       species left na = 0 while the total atom count stayed positive, so
 *       the old reader accepted it; it is now rejected at the extraction
 *       point and the warning text is asserted
 *   - LiteralZeroCountAccepted
 *     - a species with a literal zero atom count is intentionally supported
 *   - ValidFinalTripletAtExactEofAccepted
 *     - a complete final coordinate triplet ending exactly at EOF must not be
 *       misreported as a count mismatch
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
    std::vector<std::string> stru_paths;
    std::unique_ptr<Atom[]> owned_atoms;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_atoms_malformed.log";
        warn_path = testing::TempDir() + "read_atoms_malformed.warn";
        ofs_running.open(log_path);
        ofs_warning.open(warn_path);
        // a sane lattice scale: finalize_atom_positions() checks the
        // transformed coordinates, which are divided by lat0
        ucell.lat0 = 4.0;
    }

    // UnitCell::~UnitCell() deletes atoms only when set_atom_flag is true;
    // keep ownership here instead of relying on that flag
    void init_types(const int ntype)
    {
        ucell.ntype = ntype;
        owned_atoms.reset(new Atom[ntype]);
        ucell.atoms = owned_atoms.get();
        ucell.pseudo_fn.resize(ntype);
        ucell.pseudo_type.resize(ntype);
        ucell.orbital_fn.resize(ntype);
    }

    void TearDown() override
    {
        ofs_running.close();
        ofs_warning.close();
        std::remove(log_path.c_str());
        std::remove(warn_path.c_str());
        for (const std::string& path : stru_paths)
        {
            std::remove(path.c_str());
        }
    }

    std::string read_warn_file()
    {
        ofs_warning.flush();
        std::ifstream ifs(warn_path);
        return std::string(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
    }

    bool read_positions(const std::string& content)
    {
        const std::string stru_path = testing::TempDir() + "positions.stru";
        std::ofstream ofs(stru_path);
        ofs << content;
        ofs.close();
        stru_paths.push_back(stru_path);
        ofs_running << "created test STRU file: " << stru_path << std::endl;
        const std::string msg = "reading test STRU file: " + stru_path;
        std::cout << msg << std::endl;
        ofs_running << msg << std::endl;
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
    init_types(1);
    ucell.atoms[0].label = "H";
    const bool result = read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n2\n0.0 0.0 0.0\n1.0 1.0 x\n");
    EXPECT_FALSE(result);
}

TEST_F(ReadAtomsMalformedTest, NonNumericAtomCountRejected)
{
    // discriminating case: with a valid first species, the old reader left
    // the failed count at zero while the total atom count stayed positive,
    // so parsing succeeded; the new check must reject at the extraction point
    init_types(2);
    ucell.atoms[0].label = "H";
    ucell.atoms[1].label = "O";
    const bool result =
        read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n1\n0.0 0.0 0.0\nO\n0.0\nx\n");
    EXPECT_FALSE(result);
    EXPECT_THAT(read_warn_file(), testing::HasSubstr("failed to read the number of atoms"));
}

TEST_F(ReadAtomsMalformedTest, LiteralZeroCountAccepted)
{
    // a species with a literal zero atom count is intentionally supported;
    // only a failed conversion must be rejected
    init_types(2);
    ucell.atoms[0].label = "H";
    ucell.atoms[1].label = "O";
    const bool result =
        read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n1\n0.0 0.0 0.0\nO\n0.0\n0\n");
    EXPECT_TRUE(result);
    EXPECT_EQ(ucell.atoms[1].na, 0);
}

TEST_F(ReadAtomsMalformedTest, ValidFinalTripletAtExactEofAccepted)
{
    // a complete final coordinate triplet ending exactly at EOF (no trailing
    // newline) is valid input; the legacy signed-char EOF check used to
    // misreport it as "mismatch in atom number"
    init_types(1);
    ucell.atoms[0].label = "H";
    const bool result = read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n1\n0.0 0.0 0.0");
    EXPECT_TRUE(result);
    EXPECT_EQ(ucell.atoms[0].na, 1);
}

TEST_F(ReadAtomsMalformedTest, ValidPositionsAccepted)
{
    // the second atom sits at half a cell on each axis: anything on an exact
    // lattice period (e.g. 1.0 1.0 1.0) wraps onto the first atom by PBC and
    // check_tau would legitimately report overlapping atoms
    init_types(1);
    ucell.atoms[0].label = "H";
    const bool result =
        read_positions("ATOMIC_POSITIONS\nCartesian\nH\n0.0\n2\n0.0 0.0 0.0\n0.5 0.5 0.5\n");
    EXPECT_TRUE(result);
    EXPECT_EQ(ucell.atoms[0].na, 2);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].tau[1].z, 0.5);
}
