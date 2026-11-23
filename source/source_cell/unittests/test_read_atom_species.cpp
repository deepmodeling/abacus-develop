#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "source_cell/basecell.h"
#include "source_cell/read_stru.h"
#include "source_cell/unitcell.h"

/************************************************
 *  unit test of read_atom_species.cpp:
 *  named-lattice validation and malformed-input handling
 ***********************************************/

/**
 * - Tested Functions:
 *   - HexagonalMissingParamsRejected / TrigonalMissingParamsRejected
 *     - read_lattice_constant(): non-cubic named lattices require
 *       LATTICE_PARAMETERS; without it the shape variables stayed zero and
 *       produced a zero-volume or NaN cell (issue #8354)
 *   - CommentMentionDoesNotSatisfyParamsRequirement
 *     - a "# LATTICE_PARAMETERS" comment token is not a parameter section;
 *       SCAN_LINE_BEGIN matches any token in a line, so the section lookup
 *       must require the keyword to begin a line
 *   - CommentMentionDoesNotShadowRealParams
 *     - with a real LATTICE_PARAMETERS section present, the value must come
 *       from that section, not from the token after a comment mention
 *   - TrigonalMalformedParamRejected
 *     - a failed cos(alpha) extraction must not silently become the
 *       simple-cubic (alpha = 90 degree) special case
 *   - TrigonalImpossibleAngleRejected
 *     - cos(alpha) = 1.1 parses fine but yields NaN lattice vectors, which
 *       the finiteness check must reject
 *   - SmPartiallyMalformedParamsRejected / BacmPartiallyMalformedParamsRejected
 *     / TriclinicPartiallyMalformedParamsRejected
 *     - a partially malformed multi-value LATTICE_PARAMETERS line (e.g.
 *       "2 3 bad") leaves a finite nonzero-volume cell, so the determinant
 *       check alone is not enough for these branches
 *   - ValidHexagonalWithParams / ValidSmWithParams
 *     - a hexagonal lattice with c/a given still reads normally
 *   - SoRangeOverflowRejected / StRangeOverflowRejected
 *     - an out-of-range value such as 1e999 sets failbit but stores
 *       DBL_MAX, leaving a finite nonzero-volume cell that the determinant
 *       check alone would accept; every LATTICE_PARAMETERS extraction is
 *       therefore validated directly
 *   - CubicScNeedsNoParams
 *     - cubic lattice names do not need LATTICE_PARAMETERS
 *   - DegenerateLatticeVectorsRejected
 *     - explicit LATTICE_VECTORS forming a zero-volume cell are rejected
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

class ReadLatticeValidationTest : public ::testing::Test
{
  protected:
    Lattice lat;
    std::ofstream ofs_running;
    std::string log_path;
    std::vector<std::string> stru_paths;

    void SetUp() override
    {
        log_path = testing::TempDir() + "read_lattice_validation.log";
        ofs_running.open(log_path);
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

TEST_F(ReadLatticeValidationTest, HexagonalMissingParamsRejected)
{
    const std::string stru = write_stru_file("hex_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "hexagonal";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("LATTICE_PARAMETERS is required for lattice type hexagonal"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TrigonalMissingParamsRejected)
{
    const std::string stru = write_stru_file("trigonal_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "trigonal";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("LATTICE_PARAMETERS is required for lattice type trigonal"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, CommentMentionDoesNotSatisfyParamsRequirement)
{
    // a comment token "# LATTICE_PARAMETERS" is not a parameter section:
    // SCAN_LINE_BEGIN matches any token in a line, so the section lookup
    // must require the keyword to begin a line
    const std::string stru =
        write_stru_file("hex_comment_only.stru", "LATTICE_CONSTANT # LATTICE_PARAMETERS\n4.0\n");
    lat.latName = "hexagonal";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("LATTICE_PARAMETERS is required for lattice type hexagonal"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, CommentMentionDoesNotShadowRealParams)
{
    // with a real LATTICE_PARAMETERS section present, the value must come
    // from that section, not from the token after the comment mention
    const std::string stru = write_stru_file("hex_comment_plus_real.stru",
                                             "LATTICE_CONSTANT # LATTICE_PARAMETERS\n4.0\n"
                                             "LATTICE_PARAMETERS\n1.633\n");
    lat.latName = "hexagonal";
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.e33, 1.633);
    EXPECT_GT(lat.latvec.Det(), 0.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TrigonalMalformedParamRejected)
{
    // a failed cos(alpha) read must not silently become the simple-cubic
    // (alpha = 90 degree) special case
    const std::string stru =
        write_stru_file("trigonal_bad_param.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\nabc\n");
    lat.latName = "trigonal";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid cos(alpha)"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TrigonalImpossibleAngleRejected)
{
    // cos(alpha) = 1.1 parses fine but gives sqrt of a negative number, so
    // the lattice vectors become NaN; the finiteness check must reject them
    const std::string stru =
        write_stru_file("trigonal_impossible.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n1.1\n");
    lat.latName = "trigonal";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("zero-volume or non-finite cell"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, UnsupportedLatNameReportsItself)
{
    // an unsupported lattice name must report "latname not supported", not
    // be asked for LATTICE_PARAMETERS: the presence guard runs after the
    // name cascade. The legacy Warning5 fixture contains LATTICE_PARAMETERS,
    // so it passes with either guard order and cannot pin this ordering.
    const std::string stru = write_stru_file("arbitrary_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "arbitrary";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("latname not supported"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, ValidHexagonalWithParams)
{
    const std::string stru =
        write_stru_file("hex_with_param.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n1.633\n");
    lat.latName = "hexagonal";
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.e33, 1.633);
    EXPECT_GT(lat.latvec.Det(), 0.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, SmPartiallyMalformedParamsRejected)
{
    // "2 3 bad": the first two values read fine and the failed cosine read
    // leaves a finite diag(1,2,3)-like cell, which the determinant check
    // alone would accept
    const std::string stru =
        write_stru_file("sm_partial_bad.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n2 3 bad\n");
    lat.latName = "sm";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_PARAMETERS for lattice type sm"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, BacmPartiallyMalformedParamsRejected)
{
    const std::string stru =
        write_stru_file("bacm_partial_bad.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n2 3 bad\n");
    lat.latName = "bacm";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_PARAMETERS for lattice type bacm"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, TriclinicPartiallyMalformedParamsRejected)
{
    const std::string stru = write_stru_file("triclinic_partial_bad.stru",
                                             "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n0.5 0.5 bad\n");
    lat.latName = "triclinic";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_PARAMETERS for lattice type triclinic"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, ValidSmWithParams)
{
    const std::string stru =
        write_stru_file("sm_valid.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n2 3 0.5\n");
    lat.latName = "sm";
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.e21, 2.0 * 0.5);
    EXPECT_DOUBLE_EQ(lat.latvec.e33, 3.0);
    EXPECT_GT(lat.latvec.Det(), 0.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, SoRangeOverflowRejected)
{
    // "1 1e999": the out-of-range read sets failbit but stores DBL_MAX,
    // leaving det = DBL_MAX (finite and nonzero), so the determinant check
    // alone would accept this malformed input
    const std::string stru =
        write_stru_file("so_range_overflow.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n1 1e999\n");
    lat.latName = "so";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_PARAMETERS for lattice type so"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, StRangeOverflowRejected)
{
    const std::string stru =
        write_stru_file("st_range_overflow.stru", "LATTICE_CONSTANT\n4.0\nLATTICE_PARAMETERS\n1e999\n");
    lat.latName = "st";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat), ::testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("invalid LATTICE_PARAMETERS for lattice type st"));
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, CubicScNeedsNoParams)
{
    const std::string stru = write_stru_file("sc_no_param.stru", "LATTICE_CONSTANT\n4.0\n");
    lat.latName = "sc";
    std::ifstream ifa = open_stru_file(stru);
    EXPECT_TRUE(unitcell::read_lattice_constant(ifa, ofs_running, lat));
    EXPECT_DOUBLE_EQ(lat.latvec.Det(), 1.0);
    ifa.close();
}

TEST_F(ReadLatticeValidationTest, DegenerateLatticeVectorsRejected)
{
    const std::string stru = write_stru_file("degenerate.stru",
                                             "LATTICE_CONSTANT\n4.0\n"
                                             "LATTICE_VECTORS\n"
                                             "1 0 0\n"
                                             "0 1 0\n"
                                             "0 0 0\n");
    lat.latName = "user_defined_lattice";
    std::ifstream ifa = open_stru_file(stru);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(unitcell::read_lattice_constant(ifa, ofs_running, lat),
                ::testing::ExitedWithCode(1),
                "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("zero-volume"));
    ifa.close();
}

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
