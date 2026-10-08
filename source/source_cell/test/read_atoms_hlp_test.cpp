#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "source_cell/read_atoms_helper.h"
#include "source_cell/print_cell.h"
#include "source_base/vector3.h"
#include "source_base/matrix3.h"
#include "source_base/output.h"
#include <sstream>
#include <fstream>
#include <limits>

// Mock implementations for missing functions that are not in the linked sources
namespace elecstate {
    bool read_orb_file(int it, std::string& orbital_file, std::ofstream& ofs_running, Atom* atom) {
        // Mock implementation - just return true
        return true;
    }
}



// Mock Magnetism class
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

// Mock read_atom_positions function (we're testing the helpers, not the main function)
namespace unitcell {
    bool read_atom_positions(UnitCell& ucell, std::ifstream& ifpos,
                           std::ofstream& ofs_running, std::ofstream& ofs_warning) {
        // Mock implementation
        return true;
    }
}

// Test fixture for read_atoms_helper tests
class ReadAtomsHelperTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create temporary output streams
        ofs_warning.open("test_warning.log");
        ofs_running.open("test_running.log");
    }

    void TearDown() override
    {
        ofs_warning.close();
        ofs_running.close();
        // Clean up temporary files
        std::remove("test_warning.log");
        std::remove("test_running.log");
    }

    std::ofstream ofs_warning;
    std::ofstream ofs_running;
};

// Test validate_coordinate_system function
TEST_F(ReadAtomsHelperTest, ValidateCoordinateSystem_ValidInputs)
{
    EXPECT_TRUE(unitcell::validate_coordinate_system("Direct", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_angstrom", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_au", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_angstrom_center_xy", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_angstrom_center_xz", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_angstrom_center_yz", ofs_warning));
    EXPECT_TRUE(unitcell::validate_coordinate_system("Cartesian_angstrom_center_xyz", ofs_warning));
}

TEST_F(ReadAtomsHelperTest, ValidateCoordinateSystem_InvalidInputs)
{
    EXPECT_FALSE(unitcell::validate_coordinate_system("Invalid", ofs_warning));
    EXPECT_FALSE(unitcell::validate_coordinate_system("direct", ofs_warning));  // case sensitive
    EXPECT_FALSE(unitcell::validate_coordinate_system("", ofs_warning));
    EXPECT_FALSE(unitcell::validate_coordinate_system("Cartesian_angstrom_center", ofs_warning));
}

// Test calculate_lattice_center function
TEST_F(ReadAtomsHelperTest, CalculateLatticeCenterXY)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    auto center = unitcell::calculate_lattice_center(latvec, "xy");

    EXPECT_DOUBLE_EQ(center.x, 5.0);
    EXPECT_DOUBLE_EQ(center.y, 5.0);
    EXPECT_DOUBLE_EQ(center.z, 0.0);
}

TEST_F(ReadAtomsHelperTest, CalculateLatticeCenterXZ)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    auto center = unitcell::calculate_lattice_center(latvec, "xz");

    EXPECT_DOUBLE_EQ(center.x, 5.0);
    EXPECT_DOUBLE_EQ(center.y, 0.0);
    EXPECT_DOUBLE_EQ(center.z, 5.0);
}

TEST_F(ReadAtomsHelperTest, CalculateLatticeCenterYZ)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    auto center = unitcell::calculate_lattice_center(latvec, "yz");

    EXPECT_DOUBLE_EQ(center.x, 0.0);
    EXPECT_DOUBLE_EQ(center.y, 5.0);
    EXPECT_DOUBLE_EQ(center.z, 5.0);
}

TEST_F(ReadAtomsHelperTest, CalculateLatticeCenterXYZ)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    auto center = unitcell::calculate_lattice_center(latvec, "xyz");

    EXPECT_DOUBLE_EQ(center.x, 5.0);
    EXPECT_DOUBLE_EQ(center.y, 5.0);
    EXPECT_DOUBLE_EQ(center.z, 5.0);
}

TEST_F(ReadAtomsHelperTest, CalculateLatticeCenterNonCubic)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 8.0;  latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 2.0;  latvec.e22 = 6.0; latvec.e23 = 0.0;
    latvec.e31 = 1.0;  latvec.e32 = 1.0; latvec.e33 = 10.0;

    auto center = unitcell::calculate_lattice_center(latvec, "xyz");

    EXPECT_DOUBLE_EQ(center.x, (8.0 + 2.0 + 1.0) / 2.0);
    EXPECT_DOUBLE_EQ(center.y, (0.0 + 6.0 + 1.0) / 2.0);
    EXPECT_DOUBLE_EQ(center.z, (0.0 + 0.0 + 10.0) / 2.0);
}

// Test allocate_atom_properties function
TEST_F(ReadAtomsHelperTest, AllocateAtomProperties)
{
    Atom atom;
    int na = 5;
    atom.mass = 12.0;

    unitcell::allocate_atom_properties(atom, na);

    EXPECT_EQ(atom.tau.size(), na);
    EXPECT_EQ(atom.dis.size(), na);
    EXPECT_EQ(atom.taud.size(), na);
    EXPECT_EQ(atom.boundary_shift.size(), na);
    EXPECT_EQ(atom.vel.size(), na);
    EXPECT_EQ(atom.mbl.size(), na);
    EXPECT_EQ(atom.mag.size(), na);
    EXPECT_EQ(atom.angle1.size(), na);
    EXPECT_EQ(atom.angle2.size(), na);
    EXPECT_EQ(atom.m_loc_.size(), na);
    EXPECT_EQ(atom.lambda.size(), na);
    EXPECT_EQ(atom.constrain.size(), na);
    EXPECT_DOUBLE_EQ(atom.mass, 12.0);
}

// Test transform_atom_coordinates for Direct coordinates
TEST_F(ReadAtomsHelperTest, TransformAtomCoordinatesDirect)
{
    Atom atom;
    atom.tau.resize(1);
    atom.taud.resize(1);

    ModuleBase::Vector3<double> v(0.5, 0.5, 0.5);
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    double lat0 = 1.0;
    ModuleBase::Vector3<double> latcenter;

    unitcell::transform_atom_coordinates(atom, 0, "Direct", v, latvec, lat0, latcenter);

    EXPECT_DOUBLE_EQ(atom.taud[0].x, 0.5);
    EXPECT_DOUBLE_EQ(atom.taud[0].y, 0.5);
    EXPECT_DOUBLE_EQ(atom.taud[0].z, 0.5);
    EXPECT_DOUBLE_EQ(atom.tau[0].x, 5.0);
    EXPECT_DOUBLE_EQ(atom.tau[0].y, 5.0);
    EXPECT_DOUBLE_EQ(atom.tau[0].z, 5.0);
}

// Test transform_atom_coordinates for Cartesian coordinates
TEST_F(ReadAtomsHelperTest, TransformAtomCoordinatesCartesian)
{
    Atom atom;
    atom.tau.resize(1);
    atom.taud.resize(1);

    ModuleBase::Vector3<double> v(5.0, 5.0, 5.0);
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 10.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
    latvec.e21 = 0.0;  latvec.e22 = 10.0; latvec.e23 = 0.0;
    latvec.e31 = 0.0;  latvec.e32 = 0.0;  latvec.e33 = 10.0;

    double lat0 = 1.0;
    ModuleBase::Vector3<double> latcenter;

    unitcell::transform_atom_coordinates(atom, 0, "Cartesian", v, latvec, lat0, latcenter);

    EXPECT_DOUBLE_EQ(atom.tau[0].x, 5.0);
    EXPECT_DOUBLE_EQ(atom.tau[0].y, 5.0);
    EXPECT_DOUBLE_EQ(atom.tau[0].z, 5.0);
    EXPECT_DOUBLE_EQ(atom.taud[0].x, 0.5);
    EXPECT_DOUBLE_EQ(atom.taud[0].y, 0.5);
    EXPECT_DOUBLE_EQ(atom.taud[0].z, 0.5);
}

// Test process_magnetization for nspin=2
TEST_F(ReadAtomsHelperTest, ProcessMagnetizationNspin2)
{
    Atom atom;
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);

    atom.mag[0] = 2.0;
    atom.m_loc_[0].set(0, 0, 0);

    const int nspin = 2;
    const bool input_vec_mag = false;
    const bool input_angle_mag = false;
    const bool noncolin = false;
    unitcell::process_magnetization(atom, 0, 0, nspin, input_vec_mag, input_angle_mag, ofs_running, noncolin);

    // For nspin=2, only z component should be set
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].x, 0.0);
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].y, 0.0);
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].z, 2.0);
    EXPECT_DOUBLE_EQ(atom.mag[0], 2.0);
}

// Test process_magnetization for nspin=4 with vector input
TEST_F(ReadAtomsHelperTest, ProcessMagnetizationNspin4VectorInput)
{
    Atom atom;
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);

    atom.m_loc_[0].set(1.0, 1.0, 1.0);
    atom.mag[0] = sqrt(3.0);

    const int nspin = 4;
    const bool input_vec_mag = true;
    const bool input_angle_mag = false;
    const bool noncolin = true;
    unitcell::process_magnetization(atom, 0, 0, nspin, input_vec_mag, input_angle_mag, ofs_running, noncolin);

    // Angles should be calculated from vector components
    EXPECT_GT(atom.angle1[0], 0.0);
    EXPECT_GT(atom.angle2[0], 0.0);
}

// Test process_magnetization with angle input
TEST_F(ReadAtomsHelperTest, ProcessMagnetizationAngleInput)
{
    Atom atom;
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);

    atom.mag[0] = 2.0;
    atom.angle1[0] = M_PI / 2.0;  // 90 degrees
    atom.angle2[0] = 0.0;
    atom.m_loc_[0].set(0, 0, 0);

    const int nspin = 2;
    const bool input_vec_mag = false;
    const bool input_angle_mag = true;
    const bool noncolin = false;
    unitcell::process_magnetization(atom, 0, 0, nspin, input_vec_mag, input_angle_mag, ofs_running, noncolin);

    // For nspin=2, only z component is used, which should be mag[0] * cos(angle1)
    // With angle1 = PI/2, cos(PI/2) = 0
    EXPECT_NEAR(atom.m_loc_[0].z, 0.0, 1e-10);
    EXPECT_DOUBLE_EQ(atom.mag[0], atom.m_loc_[0].z);
}

// Test parse_atom_properties with movement flags
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesMovementFlags)
{
    std::string input_str = "1.0 2.0 3.0 m 1 0 1\n";
    std::istringstream iss(input_str);

    // Create a temporary file for testing
    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_EQ(mv.x, 1);
    EXPECT_EQ(mv.y, 0);
    EXPECT_EQ(mv.z, 1);

    ifpos.close();
    std::remove("test_input.tmp");
}

// Test parse_atom_properties with velocity
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesVelocity)
{
    std::string input_str = "1.0 2.0 3.0 v 0.1 0.2 0.3\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(atom.vel[0].x, 0.1);
    EXPECT_DOUBLE_EQ(atom.vel[0].y, 0.2);
    EXPECT_DOUBLE_EQ(atom.vel[0].z, 0.3);

    ifpos.close();
    std::remove("test_input.tmp");
}

// Test parse_atom_properties with scalar magnetization
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesScalarMag)
{
    std::string input_str = "1.0 2.0 3.0 mag 2.5\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(atom.mag[0], 2.5);
    EXPECT_TRUE(set_element_mag_zero);
    EXPECT_FALSE(input_vec_mag);

    ifpos.close();
    std::remove("test_input.tmp");
}

// Test parse_atom_properties with vector magnetization
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesVectorMag)
{
    std::string input_str = "1.0 2.0 3.0 mag 1.0 2.0 3.0\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].x, 1.0);
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].y, 2.0);
    EXPECT_DOUBLE_EQ(atom.m_loc_[0].z, 3.0);
    EXPECT_NEAR(atom.mag[0], sqrt(1.0 + 4.0 + 9.0), 1e-10);
    EXPECT_TRUE(input_vec_mag);
    EXPECT_TRUE(set_element_mag_zero);

    ifpos.close();
    std::remove("test_input.tmp");
}

// Test parse_atom_properties with force field (round-trip compatibility)
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesForce)
{
    std::string input_str = "1.0 2.0 3.0 m 1 1 1 f 0.5 -0.3 0.2 mag 0.8333\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_EQ(mv.x, 1);
    EXPECT_EQ(mv.y, 1);
    EXPECT_EQ(mv.z, 1);
    EXPECT_DOUBLE_EQ(atom.mag[0], 0.8333);
    EXPECT_TRUE(set_element_mag_zero);
    EXPECT_FALSE(ifpos.fail());

    ifpos.close();
    std::remove("test_input.tmp");
}

// Test parse_atom_properties with negative force values
TEST_F(ReadAtomsHelperTest, ParseAtomPropertiesNegativeForce)
{
    std::string input_str = "1.0 2.0 3.0 m 1 0 1 f -1.0 -2.0 -3.0\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "C";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(1, 1, 1);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    // Skip the position coordinates
    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_EQ(mv.x, 1);
    EXPECT_EQ(mv.y, 0);
    EXPECT_EQ(mv.z, 1);
    EXPECT_FALSE(ifpos.fail());

    ifpos.close();
    std::remove("test_input.tmp");
}

// Round-trip integration: a STRU atom line in the exact format produced by
// print_stru_file (positions + m + optional f) must be parseable by
// parse_atom_properties without leaving the stream in a fail state.
// This catches the round-trip regression reported in issue #8051.
TEST_F(ReadAtomsHelperTest, RoundTripWriterReaderForce)
{
    // Line as written by print_stru_file with has_force=true
    std::string input_str = "0.000000000000 0.000000000000 0.000000000000 m 1 1 1 f -0.123456000000 0.234567000000 -0.345678000000\n";

    std::ofstream temp_file("test_input.tmp");
    temp_file << input_str;
    temp_file.close();

    std::ifstream ifpos("test_input.tmp");

    Atom atom;
    atom.label = "Fe";
    atom.vel.resize(1);
    atom.mag.resize(1);
    atom.m_loc_.resize(1);
    atom.angle1.resize(1);
    atom.angle2.resize(1);
    atom.lambda.resize(1);
    atom.constrain.resize(1);

    ModuleBase::Vector3<int> mv(0, 0, 0);
    bool input_vec_mag = false;
    bool input_angle_mag = false;
    bool set_element_mag_zero = false;

    double x, y, z;
    ifpos >> x >> y >> z;

    bool result = unitcell::parse_atom_properties(ifpos, atom, 0, mv,
                                                  input_vec_mag, input_angle_mag,
                                                  set_element_mag_zero);

    EXPECT_TRUE(result);
    EXPECT_EQ(mv.x, 1);
    EXPECT_EQ(mv.y, 1);
    EXPECT_EQ(mv.z, 1);
    // Stream must not be in a fail state -- this was the bug in issue #8051:
    // the reader did not know "f" and consumed the force values as the next
    // keyword, leaving the stream corrupted for subsequent atoms.
    EXPECT_FALSE(ifpos.fail());
    EXPECT_TRUE(ifpos.good() || ifpos.eof());

    ifpos.close();
    std::remove("test_input.tmp");
}

// Multi-atom round-trip: write a real 3-atom spin-polarized cell with
// print_stru_file (has_force=true, nspin=2) and parse the produced file
// back through parse_atom_properties. Positions, movement flags and the
// per-atom mag value are asserted for every ia, which catches stream
// corruption that only appears when na > 1 (issue #8051).
TEST_F(ReadAtomsHelperTest, RoundTripMultiAtomMixedFields)
{
    const int nat = 3;

    // Reference data: Cartesian positions in Bohr (units of lat0),
    // movement flags and initial magnetic moments, one row per atom.
    const double pos_bohr[3][3] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.5, 0.5, 0.5}
    };
    const int mbl_int[3][3] = {
        {1, 1, 1},
        {0, 1, 0},
        {1, 0, 1}
    };
    const double mag_expected[3] = {1.5, -0.7, 2.25};

    UnitCell ucell;
    ucell.ntype = 1;
    ucell.nat = nat;
    ucell.atoms = new Atom[1];
    ucell.lat0 = 1.0;
    ucell.omega = 1.0;
    ucell.latvec.Identity();
    ucell.pseudo_fn.resize(1);
    ucell.pseudo_type.resize(1);
    ucell.pseudo_fn[0] = "Fe.upf";
    ucell.pseudo_type[0] = "uspp";
    ucell.magnet.start_mag.resize(1);
    ucell.magnet.start_mag[0] = 0.0;

    Atom& fe = ucell.atoms[0];
    fe.na = nat;
    fe.label = "Fe";
    fe.mass = 55.847;
    fe.tau.resize(nat);
    fe.mbl.resize(nat);
    fe.mag.resize(nat);
    fe.vel.resize(nat);
    fe.m_loc_.resize(nat);
    fe.angle1.resize(nat);
    fe.angle2.resize(nat);
    fe.lambda.resize(nat);
    fe.constrain.resize(nat);
    for (int ia = 0; ia < nat; ++ia)
    {
        fe.tau[ia] = ModuleBase::Vector3<double>(pos_bohr[ia][0], pos_bohr[ia][1], pos_bohr[ia][2]);
        fe.mbl[ia] = ModuleBase::Vector3<int>(mbl_int[ia][0], mbl_int[ia][1], mbl_int[ia][2]);
        fe.mag[ia] = mag_expected[ia];
    }

    // Forces in Ry/Bohr; the writer converts them to eV/Angstrom.
    ModuleBase::matrix force(nat, 3);
    for (int ia = 0; ia < nat; ++ia)
    {
        force(ia, 0) = 0.01 * (ia + 1);
        force(ia, 1) = -0.02 * (ia + 1);
        force(ia, 2) = 0.03 * (ia + 1);
    }

    const std::string filename = "test_stru_roundtrip.tmp";
    unitcell::print_stru_file(ucell, ucell.atoms, ucell.latvec, filename, "",
                              2, false, false, false, false, 0, force, true);

    std::ifstream ifpos(filename.c_str());

    // Locate the ATOMIC_POSITIONS section.
    std::string keyword;
    while (ifpos >> keyword)
    {
        if (keyword == "ATOMIC_POSITIONS")
        {
            break;
        }
    }
    ASSERT_FALSE(ifpos.fail());

    // Coordinate descriptor line.
    ifpos >> keyword;
    EXPECT_EQ(keyword, "Cartesian_angstrom");
    ifpos.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Type block header: label, default magnetism, atom count.
    std::string label;
    ifpos >> label;
    EXPECT_EQ(label, "Fe");
    ifpos.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    double header_mag = 0.0;
    ifpos >> header_mag;
    ifpos.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    int na_read = 0;
    ifpos >> na_read;
    ASSERT_EQ(na_read, nat);
    ifpos.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Parser target: every indexed buffer must cover all ia.
    Atom parse_atom;
    parse_atom.label = "Fe";
    parse_atom.vel.resize(nat);
    parse_atom.mag.resize(nat);
    parse_atom.m_loc_.resize(nat);
    parse_atom.angle1.resize(nat);
    parse_atom.angle2.resize(nat);
    parse_atom.lambda.resize(nat);
    parse_atom.constrain.resize(nat);

    const double pos_conv = ucell.lat0 * ModuleBase::BOHR_TO_A;
    for (int ia = 0; ia < nat; ++ia)
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        ifpos >> x >> y >> z;
        ASSERT_FALSE(ifpos.fail());

        ModuleBase::Vector3<int> mv(0, 0, 0);
        bool input_vec_mag = false;
        bool input_angle_mag = false;
        bool set_element_mag_zero = false;

        const bool ok = unitcell::parse_atom_properties(ifpos, parse_atom, ia, mv,
                                                        input_vec_mag, input_angle_mag,
                                                        set_element_mag_zero);
        EXPECT_TRUE(ok);

        const double expected_x = pos_bohr[ia][0] * pos_conv;
        const double expected_y = pos_bohr[ia][1] * pos_conv;
        const double expected_z = pos_bohr[ia][2] * pos_conv;
        EXPECT_NEAR(x, expected_x, 1e-8);
        EXPECT_NEAR(y, expected_y, 1e-8);
        EXPECT_NEAR(z, expected_z, 1e-8);

        EXPECT_EQ(mv.x, mbl_int[ia][0]);
        EXPECT_EQ(mv.y, mbl_int[ia][1]);
        EXPECT_EQ(mv.z, mbl_int[ia][2]);

        EXPECT_NEAR(parse_atom.mag[ia], mag_expected[ia], 1e-10);
    }

    // One read past the last consumed newline reaches end-of-file.
    ifpos.get();
    EXPECT_TRUE(ifpos.eof());
    ifpos.close();
    std::remove(filename.c_str());
    delete[] ucell.atoms;
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
