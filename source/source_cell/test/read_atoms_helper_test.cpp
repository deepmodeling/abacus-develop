#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../read_atoms_helper.h"
#include "source_base/vector3.h"
#include "source_base/matrix3.h"
#include <sstream>
#include <fstream>

// Mock implementations for missing functions that are not in the linked sources
namespace elecstate {
    bool read_orb_file(int it, std::string& orbital_file, std::ofstream& ofs_running, Atom* atom) {
        // Mock implementation - just return true
        return true;
    }
}

// Mock output class methods
void output::printM3(std::ofstream& ofs, const std::string& description, const ModuleBase::Matrix3& m) {
    // Mock implementation
}

void output::printrm(std::ofstream& ofs, const std::string& description, const ModuleBase::matrix& m, const double& limit) {
    // Mock implementation
}

// Mock InfoNonlocal class
InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}

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
    double mass = 12.0;

    unitcell::allocate_atom_properties(atom, na, mass);

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
    EXPECT_DOUBLE_EQ(atom.mass, mass);
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

    unitcell::process_magnetization(atom, 0, 0, 2, false, false, ofs_running);

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

    // Set noncolin to true to allow non-collinear magnetization
    // Note: This requires PARAM to be properly initialized

    unitcell::process_magnetization(atom, 0, 0, 4, true, false, ofs_running);

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

    // Note: For nspin=4, if noncolin is false (default), x and y components are zeroed
    // So we test with nspin=2 instead to verify the angle calculation works
    unitcell::process_magnetization(atom, 0, 0, 2, false, true, ofs_running);

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

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
