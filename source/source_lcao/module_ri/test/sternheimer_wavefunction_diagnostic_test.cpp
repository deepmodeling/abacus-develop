#include "source_lcao/module_ri/sternheimer_wavefunction_diagnostic.h"

#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <string>

namespace
{

using Diagnostic = ModuleRI::SternheimerWavefunctionDiagnostic;
using Complex = Diagnostic::Complex;

}

TEST(SternheimerWavefunctionDiagnostic, SelectorMatchesOnlyOneEquation)
{
    Diagnostic::Selector selector;
    selector.iq = 2;
    selector.ik_full = 3;
    selector.ib = 0;
    selector.ifrequency = 1;
    selector.channel = 4;

    EXPECT_TRUE(selector.matches(2, 3, 0, 1, 4));
    EXPECT_FALSE(selector.matches(1, 3, 0, 1, 4));
    EXPECT_FALSE(selector.matches(2, 2, 0, 1, 4));
    EXPECT_FALSE(selector.matches(2, 3, 1, 1, 4));
    EXPECT_FALSE(selector.matches(2, 3, 0, 2, 4));
    EXPECT_FALSE(selector.matches(2, 3, 0, 1, 5));
}

TEST(SternheimerWavefunctionDiagnostic, ParsesCompleteSelectorConfiguration)
{
    const Diagnostic::Configuration config
        = Diagnostic::parse_configuration("iq=2,ik=0,ib=0,ifreq=1,channel=4,out=wavefunction.bin");
    EXPECT_TRUE(config.selector.matches(2, 0, 0, 1, 4));
    EXPECT_EQ(config.output_filename, "wavefunction.bin");

    EXPECT_THROW(Diagnostic::parse_configuration("iq=2,ik=0,ib=0,ifreq=1,channel=4"),
                 std::invalid_argument);
    EXPECT_THROW(Diagnostic::parse_configuration("iq=2,ik=0,ib=0,ifreq=1,channel=4,out=a,iq=3"),
                 std::invalid_argument);
    EXPECT_THROW(Diagnostic::parse_configuration("iq=0,ik=0,ib=0,ifreq=1,channel=4,out=a"),
                 std::invalid_argument);
}

TEST(SternheimerWavefunctionDiagnostic, BinaryRecordRoundTripsNamedComplexVectors)
{
    Diagnostic::Record record;
    record.metadata.nx = 2;
    record.metadata.ny = 1;
    record.metadata.nz = 1;
    record.metadata.iq = 2;
    record.metadata.ik_full = 3;
    record.metadata.ib = 0;
    record.metadata.ifrequency = 1;
    record.metadata.channel = 4;
    record.metadata.lattice = {2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0};
    record.metadata.qpoint = {0.25, 0.0, 0.0};
    record.metadata.source_kpoint = {0.0, 0.0, 0.0};
    record.metadata.target_kpoint = {0.25, 0.0, 0.0};
    record.metadata.omega_ha = 0.5;
    record.metadata.omega_ry = 1.0;
    record.metadata.volume_element = 0.125;
    record.metadata.reference_eigenvalue_ry = -0.75;
    record.metadata.weighted_occupation = 0.03125;
    record.metadata.rhs_norm = 0.625;
    record.metadata.solver_relative_residual = 1.0e-9;
    record.metadata.equation_relative_residual = 2.0e-9;
    record.metadata.diagonal_branch_element = Complex(-0.125, 0.03125);
    record.vectors.push_back({"psi0", {Complex(1.0, 2.0), Complex(3.0, 4.0)}});
    record.vectors.push_back({"delta_psi", {Complex(-1.0, 0.5), Complex(0.25, -0.75)}});

    const std::string filename = "sternheimer_wavefunction_diagnostic_test.bin";
    Diagnostic::write(filename, record);
    const Diagnostic::Record reread = Diagnostic::read(filename);
    std::remove(filename.c_str());

    EXPECT_EQ(reread.metadata.nx, record.metadata.nx);
    EXPECT_EQ(reread.metadata.ny, record.metadata.ny);
    EXPECT_EQ(reread.metadata.nz, record.metadata.nz);
    EXPECT_EQ(reread.metadata.iq, record.metadata.iq);
    EXPECT_EQ(reread.metadata.ik_full, record.metadata.ik_full);
    EXPECT_EQ(reread.metadata.ib, record.metadata.ib);
    EXPECT_EQ(reread.metadata.ifrequency, record.metadata.ifrequency);
    EXPECT_EQ(reread.metadata.channel, record.metadata.channel);
    EXPECT_EQ(reread.metadata.lattice, record.metadata.lattice);
    EXPECT_EQ(reread.metadata.qpoint, record.metadata.qpoint);
    EXPECT_EQ(reread.metadata.source_kpoint, record.metadata.source_kpoint);
    EXPECT_EQ(reread.metadata.target_kpoint, record.metadata.target_kpoint);
    EXPECT_DOUBLE_EQ(reread.metadata.omega_ha, record.metadata.omega_ha);
    EXPECT_DOUBLE_EQ(reread.metadata.omega_ry, record.metadata.omega_ry);
    EXPECT_DOUBLE_EQ(reread.metadata.volume_element, record.metadata.volume_element);
    EXPECT_DOUBLE_EQ(reread.metadata.reference_eigenvalue_ry, record.metadata.reference_eigenvalue_ry);
    EXPECT_DOUBLE_EQ(reread.metadata.weighted_occupation, record.metadata.weighted_occupation);
    EXPECT_DOUBLE_EQ(reread.metadata.rhs_norm, record.metadata.rhs_norm);
    EXPECT_DOUBLE_EQ(reread.metadata.solver_relative_residual, record.metadata.solver_relative_residual);
    EXPECT_DOUBLE_EQ(reread.metadata.equation_relative_residual, record.metadata.equation_relative_residual);
    EXPECT_EQ(reread.metadata.diagonal_branch_element, record.metadata.diagonal_branch_element);
    EXPECT_EQ(reread.vectors, record.vectors);
}

TEST(SternheimerWavefunctionDiagnostic, RejectsVectorWithWrongGridSize)
{
    Diagnostic::Record record;
    record.metadata.nx = 2;
    record.metadata.ny = 2;
    record.metadata.nz = 1;
    record.vectors.push_back({"delta_psi", {Complex(1.0, 0.0)}});
    EXPECT_THROW(Diagnostic::write("unused.bin", record), std::invalid_argument);
}

TEST(SternheimerWavefunctionDiagnostic, RejectsDeclaredVectorSizeBeforeAllocation)
{
    Diagnostic::Record record;
    record.metadata.nx = 2;
    record.metadata.ny = 1;
    record.metadata.nz = 1;
    record.vectors.push_back({"psi0", {Complex(1.0, 0.0), Complex(2.0, 0.0)}});

    const std::string filename = "sternheimer_wavefunction_diagnostic_bad_size_test.bin";
    Diagnostic::write(filename, record);

    std::fstream file(filename.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    ASSERT_TRUE(file.good());
    constexpr std::streamoff vector_size_offset
        = 2 * sizeof(std::int32_t) + 8 * sizeof(std::int32_t) + 28 * sizeof(double)
          + sizeof(std::int32_t) + sizeof(std::int32_t) + 4;
    file.seekp(vector_size_offset);
    const std::int64_t invalid_size = 3;
    file.write(reinterpret_cast<const char*>(&invalid_size), sizeof(invalid_size));
    file.close();

    EXPECT_THROW(Diagnostic::read(filename), std::runtime_error);
    std::remove(filename.c_str());
}
