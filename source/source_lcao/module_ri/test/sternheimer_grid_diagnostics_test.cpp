#include "source_lcao/module_ri/sternheimer_grid_diagnostics.h"

#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using Complex = std::complex<double>;

std::string read_text_file(const std::string& path)
{
    std::ifstream input(path.c_str());
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

} // namespace

TEST(SternheimerGridDiagnostics, MeasuresResponseComponentReconstruction)
{
    const std::vector<Complex> sos = {Complex(1.0, 0.5), Complex(-0.5, 0.25)};
    const std::vector<Complex> pulay = {Complex(0.25, -0.5), Complex(0.75, 0.0)};
    const std::vector<Complex> qspace = {Complex(-0.1, 0.0), Complex(0.0, 0.5)};
    std::vector<Complex> total(sos.size());
    for (std::size_t index = 0; index != total.size(); ++index)
    {
        total[index] = sos[index] + pulay[index] + qspace[index];
    }

    EXPECT_NEAR(ModuleRI::relative_component_reconstruction_error(total, sos, pulay, qspace), 0.0, 1.0e-15);

    total[0] += Complex(1.0e-4, 0.0);
    EXPECT_GT(ModuleRI::relative_component_reconstruction_error(total, sos, pulay, qspace), 0.0);
    EXPECT_THROW(ModuleRI::relative_component_reconstruction_error(total, sos, pulay, {}), std::invalid_argument);
}

TEST(SternheimerGridDiagnostics, ValidatesProductionModeAndStableComponentFilenames)
{
    EXPECT_NO_THROW(ModuleRI::validate_grid_diagnostic_request(false, false, false, false, false));
    EXPECT_NO_THROW(ModuleRI::validate_grid_diagnostic_request(true, true, true, true, true));
    EXPECT_THROW(ModuleRI::validate_grid_diagnostic_request(true, false, true, true, true), std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_grid_diagnostic_request(true, true, false, true, true), std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_grid_diagnostic_request(true, true, true, false, true), std::invalid_argument);
    EXPECT_THROW(ModuleRI::validate_grid_diagnostic_request(true, true, true, true, false), std::invalid_argument);

    EXPECT_EQ(ModuleRI::sternheimer_component_v1_filename("sos", 1, 4, 3),
              "v1_sternheimer_component_sos_iq_1_ifreq_4_rank3.dat");
    EXPECT_THROW(ModuleRI::sternheimer_component_v1_filename("total", 1, 4, 3), std::invalid_argument);
}

TEST(SternheimerGridDiagnostics, DistributesRawKSVirtualPerturbationRows)
{
    ModuleRI::SternheimerDeltaGridFunction first_virtual;
    first_virtual.values = {Complex(1.0, 0.0), Complex(0.0, 0.0)};
    ModuleRI::SternheimerDeltaGridFunction second_virtual;
    second_virtual.values = {Complex(0.0, 0.0), Complex(0.0, 1.0)};
    const std::vector<ModuleRI::SternheimerDeltaGridFunction> virtuals = {first_virtual, second_virtual};
    const std::vector<std::vector<double>> potentials = {{2.0, 3.0}, {-1.0, 4.0}};
    const std::vector<std::vector<Complex>> occupied = {{Complex(1.0, 0.0), Complex(2.0, 0.0)}};

    auto rank0 = ModuleRI::build_local_delta_perturbation_tensor(virtuals, potentials, occupied, 0.5, 0, 2);
    const auto rank1 = ModuleRI::build_local_delta_perturbation_tensor(virtuals, potentials, occupied, 0.5, 1, 2);
    ASSERT_EQ(rank0.row_counts, (std::vector<int>{1, 0}));
    ASSERT_EQ(rank1.row_counts, (std::vector<int>{0, 1}));

    for (std::size_t index = 0; index != rank0.tensor.values.size(); ++index)
    {
        rank0.tensor.values[index] += rank1.tensor.values[index];
    }
    EXPECT_EQ(rank0.tensor.at(0, 0, 0), Complex(1.0, 0.0));
    EXPECT_EQ(rank0.tensor.at(0, 1, 0), Complex(0.0, -3.0));
    EXPECT_EQ(rank0.tensor.at(0, 0, 1), Complex(-0.5, 0.0));
    EXPECT_EQ(rank0.tensor.at(0, 1, 1), Complex(0.0, -4.0));
}

TEST(SternheimerGridDiagnostics, WritesVersionedOperatorAndPerturbationData)
{
    ModuleRI::SternheimerDeltaGridMatrices matrices;
    matrices.overlap = {Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(1.0, 0.0)};
    matrices.kinetic = {Complex(1.0, 0.0), Complex(0.1, 0.2), Complex(0.1, -0.2), Complex(2.0, 0.0)};
    matrices.local_potential = {Complex(-0.5, 0.0), Complex(0.0, 0.1), Complex(0.0, -0.1), Complex(-0.25, 0.0)};
    matrices.nonlocal = {Complex(0.2, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.3, 0.0)};
    matrices.hamiltonian.resize(4);
    for (std::size_t index = 0; index != matrices.hamiltonian.size(); ++index)
    {
        matrices.hamiltonian[index]
            = matrices.kinetic[index] + matrices.local_potential[index] + matrices.nonlocal[index];
    }

    ModuleRI::SternheimerGridDiagnosticMetadata metadata;
    metadata.nx = 2;
    metadata.ny = 1;
    metadata.nz = 1;
    metadata.spin = 1;
    metadata.occupied = 1;
    metadata.virtuals = 2;
    metadata.auxiliaries = 3;
    metadata.volume_element = 0.5;

    const std::string matrix_path = "STERNHEIMER_GRID_DIAGNOSTICS_TEST_MATRICES.dat";
    const std::string perturbation_path = "STERNHEIMER_GRID_DIAGNOSTICS_TEST_PERTURBATION.dat";
    const std::vector<Complex> occupied_virtual_overlap = {Complex(0.01, 0.02), Complex(0.03, 0.04)};

    EXPECT_NEAR(ModuleRI::relative_operator_reconstruction_error(matrices), 0.0, 1.0e-15);
    ModuleRI::write_delta_grid_matrices(matrix_path, metadata, matrices, occupied_virtual_overlap, 1.0e-12);

    const std::string matrix_text = read_text_file(matrix_path);
    EXPECT_NE(matrix_text.find("ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1"), std::string::npos);
    EXPECT_NE(matrix_text.find("grid 2 1 1"), std::string::npos);
    EXPECT_NE(matrix_text.find("spin 1"), std::string::npos);
    EXPECT_NE(matrix_text.find("energy_unit Rydberg"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix overlap"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix kinetic"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix local_potential"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix nonlocal"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix hamiltonian"), std::string::npos);
    EXPECT_NE(matrix_text.find("matrix occupied_virtual_overlap"), std::string::npos);
    EXPECT_NE(matrix_text.find("0 1"), std::string::npos);

    ModuleRI::SternheimerPerturbationTensor tensor(1, 2, 3);
    tensor.at(0, 1, 2) = Complex(0.75, -0.25);
    EXPECT_EQ(tensor.at(0, 1, 2), Complex(0.75, -0.25));
    EXPECT_THROW(tensor.at(1, 0, 0), std::out_of_range);
    ModuleRI::write_delta_perturbation_tensor(perturbation_path, metadata, tensor);

    const std::string perturbation_text = read_text_file(perturbation_path);
    EXPECT_NE(perturbation_text.find("ABACUS_STERNHEIMER_GRID_DIAGNOSTICS 1"), std::string::npos);
    EXPECT_NE(perturbation_text.find("potential_unit Rydberg"), std::string::npos);
    EXPECT_NE(perturbation_text.find("tensor perturbation occupied virtual auxiliary"), std::string::npos);
    EXPECT_NE(perturbation_text.find("0 1 2"), std::string::npos);

    ModuleRI::SternheimerDeltaGridMatrices inconsistent = matrices;
    inconsistent.hamiltonian[0] += Complex(1.0e-4, 0.0);
    EXPECT_THROW(
        ModuleRI::write_delta_grid_matrices(matrix_path, metadata, inconsistent, occupied_virtual_overlap, 1.0e-12),
        std::runtime_error);

    ModuleRI::SternheimerGridDiagnosticMetadata invalid_metadata = metadata;
    invalid_metadata.volume_element = 0.0;
    EXPECT_THROW(
        ModuleRI::write_delta_grid_matrices(matrix_path, invalid_metadata, matrices, occupied_virtual_overlap, 1.0e-12),
        std::invalid_argument);
    invalid_metadata = metadata;
    invalid_metadata.auxiliaries = 4;
    EXPECT_THROW(ModuleRI::write_delta_perturbation_tensor(perturbation_path, invalid_metadata, tensor),
                 std::invalid_argument);

    std::remove(matrix_path.c_str());
    std::remove(perturbation_path.c_str());
}
