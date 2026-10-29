#include "../snap_model.h"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>

namespace
{

TEST(SnapModel, AcceptsCompleteParameters)
{
    ModuleMD::SNAP::ModelParameters parameters;
    parameters.rcut = 4.5;
    parameters.rcutfac = 1.0;
    parameters.rmin = 0.0;
    parameters.rfac0 = 0.99363;
    parameters.twojmax = 8;
    parameters.elements.push_back("Si");
    parameters.element_weights.push_back(1.0);
    parameters.element_radii.push_back(2.25);

    std::string error;
    EXPECT_TRUE(parameters.valid(error));
    EXPECT_TRUE(error.empty());
}

TEST(SnapModel, RejectsInconsistentElementData)
{
    ModuleMD::SNAP::ModelParameters parameters;
    parameters.rcut = 4.5;
    parameters.rcutfac = 1.0;
    parameters.rmin = 0.0;
    parameters.rfac0 = 0.99363;
    parameters.elements.push_back("Si");
    parameters.element_radii.push_back(2.25);

    std::string error;
    EXPECT_FALSE(parameters.valid(error));
    EXPECT_FALSE(error.empty());
}

TEST(SnapModel, RejectsInvalidRadialRange)
{
    ModuleMD::SNAP::ModelParameters parameters;
    parameters.rcut = 3.0;
    parameters.rcutfac = 1.0;
    parameters.rmin = 3.0;
    parameters.rfac0 = 1.0;
    parameters.elements.push_back("Si");
    parameters.element_weights.push_back(1.0);
    parameters.element_radii.push_back(2.25);

    std::string error;
    EXPECT_FALSE(parameters.valid(error));
    EXPECT_FALSE(error.empty());
}

TEST(SnapModel, LoadsLammpsStyleCoefficientAndParameterFiles)
{
    const std::string coefficient_file = testing::TempDir() + "snap_coeff_test.txt";
    const std::string parameter_file = testing::TempDir() + "snap_param_test.txt";
    {
        std::ofstream output(coefficient_file.c_str());
        output << "1 3\n";
        output << "Si 2.25 1.0\n";
        output << "-0.5\n1.25\n2.5\n";
    }
    {
        std::ofstream output(parameter_file.c_str());
        output << "rcutfac 1.2\n";
        output << "twojmax 4\n";
        output << "rfac0 0.99363\n";
        output << "rmin0 0.0\n";
        output << "bnormflag 0\n";
        output << "chemflag 0\n";
    }

    ModuleMD::SNAP::ModelData model;
    std::string error;
    EXPECT_TRUE(model.load(coefficient_file, parameter_file, std::vector<std::string>(1, "Si"), error));
    EXPECT_TRUE(error.empty());
    ASSERT_EQ(model.coefficient_count(), static_cast<std::size_t>(3));
    EXPECT_DOUBLE_EQ(model.parameters.element_weights[0], 1.0);
    EXPECT_DOUBLE_EQ(model.parameters.element_radii[0], 2.25);
    EXPECT_DOUBLE_EQ(model.parameters.rcut, 5.4);
    EXPECT_DOUBLE_EQ(model.coefficients[0][1], 1.25);

    std::remove(coefficient_file.c_str());
    std::remove(parameter_file.c_str());
}

TEST(SnapModel, LoadsMultipleElementBlocks)
{
    const std::string coefficient_file = testing::TempDir() + "snap_coeff_multi_test.txt";
    const std::string parameter_file = testing::TempDir() + "snap_param_multi_test.txt";
    {
        std::ofstream output(coefficient_file.c_str());
        output << "2 2\n";
        output << "Si 2.25 1.0\n";
        output << "-0.5\n1.25\n";
        output << "C 1.7 0.8\n";
        output << "-0.25\n2.5\n";
    }
    {
        std::ofstream output(parameter_file.c_str());
        output << "rcutfac 1.0\n";
        output << "twojmax 0\n";
        output << "rfac0 0.99363\n";
        output << "rmin0 0.0\n";
        output << "chemflag 0\n";
    }

    ModuleMD::SNAP::ModelData model;
    std::string error;
    const std::vector<std::string> elements = {"Si", "C"};
    EXPECT_TRUE(model.load(coefficient_file, parameter_file, elements, error));
    EXPECT_TRUE(error.empty());
    ASSERT_EQ(model.coefficients.size(), static_cast<std::size_t>(2));
    EXPECT_DOUBLE_EQ(model.parameters.element_weights[1], 0.8);
    EXPECT_DOUBLE_EQ(model.parameters.element_radii[1], 1.7);
    EXPECT_DOUBLE_EQ(model.coefficients[1][1], 2.5);

    std::remove(coefficient_file.c_str());
    std::remove(parameter_file.c_str());
}

TEST(SnapModel, LoadsChemicalCoefficientBlocks)
{
    const std::string coefficient_file = testing::TempDir() + "snap_coeff_chemical_test.txt";
    const std::string parameter_file = testing::TempDir() + "snap_param_chemical_test.txt";
    {
        std::ofstream output(coefficient_file.c_str());
        output << "2 9\n";
        output << "Si 2.25 1.0\n";
        for (int i = 0; i < 9; ++i)
        {
            output << (i == 1 ? 1.0 : 0.0) << "\n";
        }
        output << "C 1.7 0.8\n";
        for (int i = 0; i < 9; ++i)
        {
            output << (i == 1 ? 2.0 : 0.0) << "\n";
        }
    }
    {
        std::ofstream output(parameter_file.c_str());
        output << "rcutfac 1.0\n";
        output << "twojmax 0\n";
        output << "rfac0 0.99363\n";
        output << "rmin0 0.0\n";
        output << "chemflag 1\n";
        output << "bnormflag 1\n";
        output << "wselfallflag 1\n";
    }

    ModuleMD::SNAP::ModelData model;
    std::string error;
    const std::vector<std::string> elements = {"Si", "C"};
    EXPECT_TRUE(model.load(coefficient_file, parameter_file, elements, error));
    EXPECT_TRUE(error.empty());
    ASSERT_EQ(model.coefficients.size(), static_cast<std::size_t>(2));
    EXPECT_TRUE(model.parameters.chem_flag);
    EXPECT_EQ(model.coefficients[0].size(), static_cast<std::size_t>(9));
    EXPECT_DOUBLE_EQ(model.coefficients[1][1], 2.0);

    std::remove(coefficient_file.c_str());
    std::remove(parameter_file.c_str());
}

} // namespace
