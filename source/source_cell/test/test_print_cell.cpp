#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_cell/cal_ux.h"
#include "source_cell/read_orb.h"
#include "source_cell/read_pp_ucell.h"
#include "source_cell/read_stru.h"
#include "source_cell/cell_tools.h"
#include "source_cell/print_cell.h"
#include "memory"
#include "source_base/global_variable.h"
#include "source_base/mathzone.h"
#include "prepare_unitcell.h"
#include "source_cell/update_cell.h"
#include <fstream>
#include <cstdio>
#include <streambuf>
#include <valarray>
#include <vector>

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

class PrintCellTest : public testing::Test
{
protected:
    std::unique_ptr<UnitCell> ucell;
    void SetUp() override
    {
        // nothing to do here, each test sets up its own ucell
    }
    void TearDown() override
    {
        // cleanup generated files
    }
};

TEST_F(PrintCellTest, PrintSTRU_nspin1)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin1";

    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 1, false, false, false, false, false, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ATOMIC_SPECIES"));
    EXPECT_THAT(str, testing::HasSubstr("C  12.0000 C.upf upf201"));
    EXPECT_THAT(str, testing::HasSubstr("H   1.0000 H.upf upf201"));
    EXPECT_THAT(str, testing::HasSubstr("LATTICE_CONSTANT"));
    // lat0 is always printed as one Angstrom in Bohr
    EXPECT_THAT(str, testing::HasSubstr(
        "1.8897268778 # in Bohr (= 1 Angstrom); lattice vectors below are in Angstrom"));
    EXPECT_THAT(str, testing::HasSubstr("LATTICE_VECTORS # in Angstrom"));
    EXPECT_THAT(str, testing::HasSubstr("ATOMIC_POSITIONS"));
    // Cartesian (direct=false) is always emitted as Cartesian_angstrom
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom # positions in Angstrom"));
    EXPECT_THAT(str, testing::HasSubstr("C #label"));
    EXPECT_THAT(str, testing::HasSubstr("0.0000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str, testing::HasSubstr("1 #number of atoms"));
    // single-space separated coordinates
    EXPECT_THAT(str, testing::HasSubstr("0.9999996019 0.9999996019 0.9999996019 m 1 1 1"));
    EXPECT_THAT(str, testing::HasSubstr("H #label"));
    EXPECT_THAT(str, testing::HasSubstr("0.0000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str, testing::HasSubstr("2 #number of atoms"));
    EXPECT_THAT(str, testing::HasSubstr("1.4999994028 1.4999994028 1.4999994028 m 0 0 0"));
    EXPECT_THAT(str, testing::HasSubstr("0.4999998009 0.4999998009 0.4999998009 m 0 0 1"));
    // No force output when force matrix is empty
    EXPECT_THAT(str, testing::Not(testing::HasSubstr(" f ")));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_nspin2_no_force)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin2";

    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 2, true, true, false, false, false, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("ATOMIC_POSITIONS"));
    // direct=true without forces keeps fractional coordinates
    EXPECT_THAT(str, testing::HasSubstr("Direct\n"));
    EXPECT_THAT(str, testing::HasSubstr("C #label"));
    EXPECT_THAT(str, testing::HasSubstr("0.0000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1000000000 0.1000000000 0.1000000000 m 1 1 1 v 0.1000000000 0.1000000000 0.1000000000 mag 0.0000"));
    EXPECT_THAT(str, testing::HasSubstr("H #label"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1500000000 0.1500000000 0.1500000000 m 0 0 0 v 0.1000000000 0.1000000000 0.1000000000 mag 0.0000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.0500000000 0.0500000000 0.0500000000 m 0 0 1 v 0.1000000000 0.1000000000 0.1000000000 mag 0.0000"));
    // No force output when force matrix is empty
    EXPECT_THAT(str, testing::Not(testing::HasSubstr(" f ")));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_nspin2_with_mag)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin2_mag";

    ucell->atoms[0].mag[0] = 1.5;   // C
    ucell->atoms[1].mag[0] = -0.5;  // H1
    ucell->atoms[1].mag[1] = 2.0;   // H2
    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 2, true, false, false, false, false, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("C #label"));
    EXPECT_THAT(str, testing::HasSubstr("1.5000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1000000000 0.1000000000 0.1000000000 m 1 1 1 mag 1.5000"));
    EXPECT_THAT(str, testing::HasSubstr("H #label"));
    EXPECT_THAT(str, testing::HasSubstr("-0.5000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1500000000 0.1500000000 0.1500000000 m 0 0 0 mag -0.5000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.0500000000 0.0500000000 0.0500000000 m 0 0 1 mag 2.0000"));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_nspin2_mulliken)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin2_mulliken";

    ucell->descriptor_file = "__unittest_numerical_descriptor__";
    ucell->orbital_fn[0] = "__unittest_orbital_fn_0__";
    ucell->orbital_fn[1] = "__unittest_orbital_fn_1__";
    ucell->atom_mulliken = {{-1, 0.5}, {-1, 0.4}, {-1, 0.3}};
    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 2, true, false, true, true, true, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("NUMERICAL_ORBITAL"));
    EXPECT_THAT(str, testing::HasSubstr("__unittest_orbital_fn_0__"));
    EXPECT_THAT(str, testing::HasSubstr("__unittest_orbital_fn_1__"));
    EXPECT_THAT(str, testing::HasSubstr("NUMERICAL_DESCRIPTOR"));
    EXPECT_THAT(str, testing::HasSubstr("__unittest_numerical_descriptor__"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1000000000 0.1000000000 0.1000000000 m 1 1 1 mag 0.5000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1500000000 0.1500000000 0.1500000000 m 0 0 0 mag 0.4000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.0500000000 0.0500000000 0.0500000000 m 0 0 1 mag 0.3000"));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_nspin4_initial_mag)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin4";

    ucell->atoms[0].m_loc_[0].set(1.0, 0.0, 0.0);   // C
    ucell->atoms[1].m_loc_[0].set(0.0, 1.0, 0.0);   // H1
    ucell->atoms[1].m_loc_[1].set(0.0, 0.0, 1.0);   // H2
    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 4, true, false, false, false, false, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str, testing::HasSubstr("C #label"));
    EXPECT_THAT(str, testing::HasSubstr("1.0000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1000000000 0.1000000000 0.1000000000 m 1 1 1 mag 1.0000 0.0000 0.0000"));
    EXPECT_THAT(str, testing::HasSubstr("H #label"));
    EXPECT_THAT(str, testing::HasSubstr("1.0000 #magnetism (default, overridden by per-atom mag below)"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1500000000 0.1500000000 0.1500000000 m 0 0 0 mag 0.0000 1.0000 0.0000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.0500000000 0.0500000000 0.0500000000 m 0 0 1 mag 0.0000 0.0000 1.0000"));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_nspin4_mulliken)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_nspin4_mulliken";

    ucell->atom_mulliken = {{-1, 0.5, 0.1, 0.2}, {-1, 0.4, 0.3, 0.4}, {-1, 0.3, 0.5, 0.6}};
    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 4, true, false, true, false, false, 0, ModuleBase::matrix());
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    EXPECT_THAT(str,
                testing::HasSubstr("0.1000000000 0.1000000000 0.1000000000 m 1 1 1 mag 0.5000 0.1000 0.2000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.1500000000 0.1500000000 0.1500000000 m 0 0 0 mag 0.4000 0.3000 0.4000"));
    EXPECT_THAT(str,
                testing::HasSubstr("0.0500000000 0.0500000000 0.0500000000 m 0 0 1 mag 0.3000 0.5000 0.6000"));
    ifs.close();
    remove(fn.c_str());
}

TEST_F(PrintCellTest, PrintSTRU_with_force)
{
    UcellTestPrepare utp = UcellTestLib["C1H2-Index"];
    ucell = utp.SetUcellInfo();
    std::string fn = "C1H2_STRU_force";

    // Create force matrix: nat=3, nc=3 (internal unit Ry/Bohr)
    ModuleBase::matrix force(3, 3);
    force(0, 0) = 0.1; force(0, 1) = 0.2; force(0, 2) = 0.3;    // C
    force(1, 0) = -0.1; force(1, 1) = -0.2; force(1, 2) = -0.3; // H1
    force(2, 0) = 0.05; force(2, 1) = 0.15; force(2, 2) = -0.25; // H2

    unitcell::print_stru_file(*ucell, ucell->atoms, ucell->latvec,
                              fn, "", 2, true, false, false, false, false, 0, force);
    std::ifstream ifs;
    ifs.open(fn);
    std::string str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // With forces fractional coordinates are replaced by Cartesian_angstrom
    EXPECT_THAT(str, testing::HasSubstr("Cartesian_angstrom # positions in Angstrom, forces in eV/Angstrom"));

    // Internal Cartesian tau is in units of lat0 (Bohr); expected Angstrom = tau * lat0 * BOHR_TO_A
    const double pos_conv = ucell->lat0 * ModuleBase::BOHR_TO_A;
    // Internal force is Ry/Bohr; expected eV/Angstrom = f * Ry_to_eV / BOHR_TO_A
    const double force_conv = ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A;

    auto fmt3 = [](double a, double b, double c) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "%.10f %.10f %.10f", a, b, c);
        return std::string(buf);
    };
    auto fmtf = [](double a, double b, double c) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), " f %.6f %.6f %.6f", a, b, c);
        return std::string(buf);
    };

    // C (iat=0): m 1 1 1
    EXPECT_THAT(str, testing::HasSubstr(
        fmt3(ucell->atoms[0].tau[0].x * pos_conv,
             ucell->atoms[0].tau[0].y * pos_conv,
             ucell->atoms[0].tau[0].z * pos_conv) + " m 1 1 1"));
    EXPECT_THAT(str, testing::HasSubstr(
        fmtf(force(0, 0) * force_conv, force(0, 1) * force_conv, force(0, 2) * force_conv)));

    // H1 (iat=1): m 0 0 0
    EXPECT_THAT(str, testing::HasSubstr(
        fmt3(ucell->atoms[1].tau[0].x * pos_conv,
             ucell->atoms[1].tau[0].y * pos_conv,
             ucell->atoms[1].tau[0].z * pos_conv) + " m 0 0 0"));
    EXPECT_THAT(str, testing::HasSubstr(
        fmtf(force(1, 0) * force_conv, force(1, 1) * force_conv, force(1, 2) * force_conv)));

    // H2 (iat=2): m 0 0 1
    EXPECT_THAT(str, testing::HasSubstr(
        fmt3(ucell->atoms[1].tau[1].x * pos_conv,
             ucell->atoms[1].tau[1].y * pos_conv,
             ucell->atoms[1].tau[1].z * pos_conv) + " m 0 0 1"));
    EXPECT_THAT(str, testing::HasSubstr(
        fmtf(force(2, 0) * force_conv, force(2, 1) * force_conv, force(2, 2) * force_conv)));

    ifs.close();
    remove(fn.c_str());
}
