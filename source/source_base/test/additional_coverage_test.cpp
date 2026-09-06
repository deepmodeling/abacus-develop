#include "source_base/clebsch_gordan_coeff.h"
#include "source_base/cubic_spline.h"
#include "source_base/global_function.h"
#include "source_base/math_integral.h"
#include "source_base/math_lebedev_laikov.h"
#include "source_base/mathzone_add1.h"
#include "source_base/matrix.h"
#include "source_base/module_device/device.h"
#include "source_base/module_mixing/mixing_data.h"
#include "source_base/module_out/file_reader.h"
#include "source_base/output.h"
#include "source_base/tool_quit.h"
#include "source_base/vector3.h"
#include "source_base/ylm.h"

#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

TEST(SourceBaseAdditionalCoverage, ClebschGordanLifecycle)
{
    ModuleBase::Clebsch_Gordan clebsch_gordan;
}

TEST(SourceBaseAdditionalCoverage, MatrixMove)
{
    ModuleBase::matrix source_matrix(2, 2, true);
    source_matrix(0, 0) = 3.0;
    ModuleBase::matrix moved_matrix(std::move(source_matrix));
    EXPECT_EQ(moved_matrix.nr, 2);
    EXPECT_DOUBLE_EQ(moved_matrix(0, 0), 3.0);
}

TEST(SourceBaseAdditionalCoverage, VectorConstructors)
{
    const int x = 1;
    const int y = 2;
    const int z = 3;
    ModuleBase::Vector3<int> integer_vector(x, y, z);
    EXPECT_EQ(integer_vector.z, 3);

    ModuleBase::Vector3<double> source_vector(1.0, 2.0, 3.0);
    ModuleBase::Vector3<double> moved_vector(std::move(source_vector));
    EXPECT_DOUBLE_EQ(moved_vector.y, 2.0);
}

TEST(SourceBaseAdditionalCoverage, MixingDataLifecycle)
{
    Base_Mixing::Mixing_Data mixing_data(2, 3, sizeof(double));
    EXPECT_NE(mixing_data.data, nullptr);
    EXPECT_EQ(mixing_data.ndim_tot, 2);
    EXPECT_EQ(mixing_data.length, 3);
}

TEST(SourceBaseAdditionalCoverage, NumericalWrappers)
{
    const int count = 3;
    double points[count] = {};
    double weights[count] = {};
    ModuleBase::Integral::Gauss_Legendre_grid_and_weight(count, points, weights);
    EXPECT_NEAR(points[0], -std::sqrt(3.0 / 5.0), 1.0e-12);
    EXPECT_NEAR(points[1], 0.0, 1.0e-12);
    EXPECT_NEAR(weights[0] + weights[1] + weights[2], 2.0, 1.0e-12);

    double scaled_points[count] = {};
    double scaled_weights[count] = {};
    ModuleBase::Integral::Gauss_Legendre_grid_and_weight(0.0, 2.0, count, scaled_points, scaled_weights);
    EXPECT_NEAR(scaled_points[1], 1.0, 1.0e-12);
    EXPECT_NEAR(scaled_weights[0] + scaled_weights[1] + scaled_weights[2], 2.0, 1.0e-12);

    const std::complex<float> left[2] = {{1.0F, 2.0F}, {3.0F, 4.0F}};
    const std::complex<float> right[2] = {{2.0F, 1.0F}, {4.0F, 3.0F}};
    EXPECT_FLOAT_EQ(ModuleBase::GlobalFunc::ddot_real(2, left, right, false), 28.0F);

    const double knots[3] = {0.0, 1.0, 2.0};
    ModuleBase::CubicSpline spline(3, knots);
    EXPECT_DOUBLE_EQ(spline.xmin(), 0.0);
    EXPECT_DOUBLE_EQ(spline.xmax(), 2.0);

    const double radial_values[3] = {0.0, 1.0, 4.0};
    double derivative[3] = {};
    ModuleBase::Mathzone_Add1::Uni_Deriv_Phi(radial_values, 3, 1.0, 1, derivative);
    EXPECT_TRUE(std::isfinite(derivative[1]));
}

TEST(SourceBaseAdditionalCoverage, LegacySphericalHarmonics)
{
    const int lmax = 2;
    const ModuleBase::Vector3<double> direction(1.0, 0.0, 0.0);
    double values[4] = {};
    double gradients[4][3] = {};
    ModuleBase::Ylm::get_ylm_real(lmax, direction, values, gradients);
    EXPECT_TRUE(std::isfinite(values[0]));
    EXPECT_TRUE(std::isfinite(gradients[1][0]));

    double solid_values[4] = {};
    ModuleBase::Ylm::rlylm(lmax, 1.0, 0.0, 0.0, solid_values);
    EXPECT_TRUE(std::isfinite(solid_values[0]));
}

TEST(SourceBaseAdditionalCoverage, OutputAndConfigurationHelpers)
{
    const std::string output_file = "source_base_additional_coverage.log";
    std::ofstream output_stream(output_file.c_str());
    const ModuleBase::Matrix3 matrix(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    output::printM3(output_stream, "matrix", matrix);
    base_device::information::record_device_memory<base_device::DEVICE_CPU>(nullptr, output_stream, "cpu", 0);
    output_stream.close();

    std::ifstream input_stream(output_file.c_str());
    std::string line;
    std::getline(input_stream, line);
    EXPECT_NE(line.find("matrix"), std::string::npos);
    input_stream.close();
    std::remove(output_file.c_str());

    ModuleBase::set_quit_out_dir("coverage-output/");
    EXPECT_EQ(ModuleBase::get_global_out_dir(), "coverage-output/");
    ModuleBase::set_quit_calculation("unit-test");
    ModuleBase::CHECK_WARNING_QUIT(false, "coverage", "no error");
    ModuleBase::GlobalFunc::NOTE("covered no-op");
}

TEST(SourceBaseAdditionalCoverage, UnitCellReader)
{
    const std::string unit_cell_file = "source_base_unit_cell_coverage.txt";
    std::ofstream unit_cell_output(unit_cell_file.c_str());
    unit_cell_output << "lattice name\n";
    unit_cell_output << "1.0\n";
    unit_cell_output << "1 0 0\n";
    unit_cell_output << "0 1 0\n";
    unit_cell_output << "0 0 1\n";
    unit_cell_output << "H He\n";
    unit_cell_output << "1 1\n";
    unit_cell_output << "Direct\n";
    unit_cell_output << "0 0 0\n";
    unit_cell_output << "0.5 0.5 0.5\n";
    unit_cell_output << "after unit cell\n";
    unit_cell_output.close();

    ModuleIO::FileReader unit_cell_reader(unit_cell_file);
    unit_cell_reader.read_ucell();
    unit_cell_reader.readLine();
    EXPECT_EQ(unit_cell_reader.ss.str(), "after unit cell");
    std::remove(unit_cell_file.c_str());
}

TEST(SourceBaseAdditionalCoverage, FileScanningAndLebedevOutput)
{
    const std::string scan_file = "source_base_scan_coverage.txt";
    std::ofstream scan_output(scan_file.c_str());
    scan_output << "# target in a comment\n";
    scan_output << "prefix target suffix\n";
    scan_output.close();

    std::ifstream scan_input(scan_file.c_str());
    EXPECT_TRUE(ModuleBase::GlobalFunc::SCAN_LINE_BEGIN(scan_input, "target", true, false));
    scan_input.close();
    std::remove(scan_file.c_str());

    ModuleBase::Lebedev_laikov_grid grid(6);
    grid.generate_grid_points();
    grid.print_grid_and_weight("source_base_lebedev_coverage");
    std::ifstream lebedev_output("source_base_lebedev_coverage_degree6");
    EXPECT_TRUE(lebedev_output.good());
    lebedev_output.close();
    std::remove("source_base_lebedev_coverage_degree6");
}
