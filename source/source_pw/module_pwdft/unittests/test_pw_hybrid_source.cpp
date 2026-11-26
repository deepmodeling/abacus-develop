#include "source_pw/module_pwdft/pw_hybrid_source.h"
#include "source_base/matrix.h"

#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace
{
std::vector<pw::ExxSourceHeader> source_headers(int nspin)
{
    const int nks = nspin == 2 ? 4 : 2;
    std::vector<pw::ExxSourceHeader> headers(nks);
    for (int iq = 0; iq < nks; ++iq)
    {
        auto& header = headers[iq];
        header.ik = iq + 1;
        header.nks = nks;
        header.npw = 100;
        header.nbands = 3;
        header.kvec_c[0] = 0.123456789123456 * (iq % 2);
        header.weight = nspin == 1 ? 1.0 : 0.5;
    }
    return headers;
}

std::string eig_occ(int nspin)
{
    const std::vector<pw::ExxSourceHeader> headers = source_headers(nspin);
    std::ostringstream out;
    out << "# ABACUS version: test\n# Written at 2026-01-01 00:00:00\n# Ionic_Step 1\n"
        << " Electronic state energy (eV) and occupations\n Spin number " << nspin << '\n';
    for (int iq = 0; iq < headers.size(); ++iq)
    {
        const auto& header = headers[iq];
        out << " spin=" << iq / 2 + 1 << " k-point=" << iq % 2 + 1
            << "/2 Cartesian=" << std::setprecision(8) << header.kvec_c[0]
            << " 0 0 (100 plane wave)\n";
        out << "1 -1 " << header.weight << "\n2 0 " << header.weight * 0.125 << "\n3 1 0\n\n";
    }
    return out.str();
}
}

TEST(PWHybridSource, ReadsExistingWeightedOccupationsForBothSpins)
{
    for (int nspin = 1; nspin <= 2; ++nspin)
    {
        const std::vector<pw::ExxSourceHeader> headers = source_headers(nspin);
        std::istringstream in(eig_occ(nspin));
        ModuleBase::matrix weights;
        pw::read_exx_source_occupations(in, headers, nspin, weights);
        EXPECT_EQ(weights.nr, headers.size());
        EXPECT_EQ(weights.nc, 3);
        EXPECT_DOUBLE_EQ(weights(1, 1), headers[1].weight * 0.125);
    }
}

TEST(PWHybridSource, RejectsInvalidExistingOccupationFiles)
{
    const std::vector<pw::ExxSourceHeader> headers = source_headers(1);
    ModuleBase::matrix weights;
    const std::string saved = eig_occ(1);
    std::istringstream truncated(saved.substr(0, saved.size() - 6));
    EXPECT_THROW(pw::read_exx_source_occupations(truncated, headers, 1, weights), std::runtime_error);
    std::istringstream trailing(saved + saved);
    EXPECT_THROW(pw::read_exx_source_occupations(trailing, headers, 1, weights), std::runtime_error);
    std::istringstream wrong_spin(saved);
    EXPECT_THROW(pw::read_exx_source_occupations(wrong_spin, headers, 2, weights), std::runtime_error);
    std::vector<pw::ExxSourceHeader> changed = headers;
    changed[1].kvec_c[0] += 0.1;
    std::istringstream wrong_coordinates(saved);
    EXPECT_THROW(pw::read_exx_source_occupations(wrong_coordinates, changed, 1, weights), std::runtime_error);
    std::string invalid = saved;
    invalid.replace(invalid.find("0.125"), 5, "-0.25");
    std::istringstream negative(invalid);
    try
    {
        // A negative occupation is what a Methfessel-Paxton source can produce; the
        // diagnostic must name the unsupported smearing rather than report invalid data.
        pw::read_exx_source_occupations(negative, headers, 1, weights);
        FAIL() << "Expected std::runtime_error";
    }
    catch (const std::runtime_error& error)
    {
        EXPECT_NE(std::string(error.what()).find("Methfessel-Paxton"), std::string::npos);
        EXPECT_NE(std::string(error.what()).find("smearing_method fixed"), std::string::npos);
    }
}

TEST(PWHybridSource, RejectsDifferentMeshesInSpinChannels)
{
    std::vector<pw::ExxSourceHeader> headers = source_headers(2);
    headers[2].kvec_c[0] = 0.2;
    std::string saved = eig_occ(2);
    const std::string original = "spin=2 k-point=1/2 Cartesian=0";
    const std::size_t location = saved.find(original);
    ASSERT_NE(location, std::string::npos);
    saved.replace(location, original.size(), "spin=2 k-point=1/2 Cartesian=0.2");
    std::istringstream in(saved);
    ModuleBase::matrix weights;
    EXPECT_THROW(pw::read_exx_source_occupations(in, headers, 2, weights), std::runtime_error);
}

TEST(PWHybridSource, ReadsExistingBinaryWavefunctionHeader)
{
    // Native binary record layout from write_wfc_pw; no new wavefunction format.
    std::ostringstream out(std::ios::binary);
    const int record = 72;
    const int ik = 1;
    const int nks = 2;
    const int npw = 1;
    const int nbands = 3;
    const double coordinate = 0.123456789123456;
    const double weight = 1;
    const double cutoff = 10;
    const double lattice = 5;
    const double tpiba = 1.256637061435917;
    out.write(reinterpret_cast<const char*>(&record), sizeof(record));
    out.write(reinterpret_cast<const char*>(&ik), sizeof(ik));
    out.write(reinterpret_cast<const char*>(&nks), sizeof(nks));
    for (int axis = 0; axis < 3; ++axis)
    {
        out.write(reinterpret_cast<const char*>(&coordinate), sizeof(coordinate));
    }
    out.write(reinterpret_cast<const char*>(&weight), sizeof(weight));
    out.write(reinterpret_cast<const char*>(&npw), sizeof(npw));
    out.write(reinterpret_cast<const char*>(&nbands), sizeof(nbands));
    out.write(reinterpret_cast<const char*>(&cutoff), sizeof(cutoff));
    out.write(reinterpret_cast<const char*>(&lattice), sizeof(lattice));
    out.write(reinterpret_cast<const char*>(&tpiba), sizeof(tpiba));
    out.write(reinterpret_cast<const char*>(&record), sizeof(record));
    std::string saved = out.str();
    saved.resize(160 + 8 + 12 * npw + nbands * (8 + 16 * npw), '\0');
    std::istringstream in(saved, std::ios::binary);
    const pw::ExxSourceHeader header = pw::read_exx_source_header(in);
    EXPECT_EQ(header.nbands, 3);
    EXPECT_DOUBLE_EQ(header.kvec_c[0], coordinate);
    EXPECT_DOUBLE_EQ(header.weight, weight);
    std::istringstream truncated(saved.substr(0, saved.size() - 1), std::ios::binary);
    EXPECT_THROW(pw::read_exx_source_header(truncated), std::runtime_error);
    saved[0] = 0;
    std::istringstream invalid(saved, std::ios::binary);
    EXPECT_THROW(pw::read_exx_source_header(invalid), std::runtime_error);
}

namespace
{
/**
 * @brief Create an empty placeholder file for a resolver test.
 */
void touch(const std::string& path)
{
    std::ofstream out(path.c_str(), std::ios::binary);
    out << "source";
}

/**
 * @brief Resolve a name inside the gtest temporary directory.
 */
std::string temp_path(const std::string& name)
{
    return testing::TempDir() + name;
}
} // namespace

TEST(PWHybridSource, ResolverPrefersFixedName)
{
    const std::string fixed = temp_path("resolver_fixed.txt");
    const std::string numbered = temp_path("resolver_fixedg2.txt");
    touch(fixed);
    touch(numbered);
    const std::vector<std::string> directories(1, testing::TempDir());
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_fixed.txt", ".txt", false), fixed);
    std::remove(fixed.c_str());
    std::remove(numbered.c_str());
}

TEST(PWHybridSource, ResolverSelectsLatestIonicStep)
{
    const std::vector<std::string> directories(1, testing::TempDir());
    const std::string g1 = temp_path("resolver_stepg1.txt");
    const std::string g3 = temp_path("resolver_stepg3.txt");
    touch(g3);
    touch(g1);
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_step.txt", ".txt", false), g3);
    std::remove(g1.c_str());
    std::remove(g3.c_str());
}

TEST(PWHybridSource, ResolverSelectsLatestElectronicStepWithinLatestIonicStep)
{
    const std::vector<std::string> directories(1, testing::TempDir());
    const std::string g1e7 = temp_path("resolver_wfg1e7_pw.dat");
    const std::string g2e3 = temp_path("resolver_wfg2e3_pw.dat");
    const std::string g2e5 = temp_path("resolver_wfg2e5_pw.dat");
    touch(g1e7);
    touch(g2e3);
    touch(g2e5);
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_wf_pw.dat", "_pw.dat", true), g2e5);
    std::remove(g1e7.c_str());
    std::remove(g2e3.c_str());
    std::remove(g2e5.c_str());
}

TEST(PWHybridSource, ResolverReportsSnapshotTokens)
{
    const std::vector<std::string> directories(1, testing::TempDir());
    // Fixed (step-less) name reports -1/-1 so callers can lock any snapshot.
    const std::string fixed = temp_path("resolver_snap.txt");
    touch(fixed);
    int ionic_step = 0;
    int electronic_step = 0;
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_snap.txt", ".txt", false,
                                      &ionic_step, &electronic_step), fixed);
    EXPECT_EQ(ionic_step, -1);
    EXPECT_EQ(electronic_step, -1);
    std::remove(fixed.c_str());
    // g{step}e{iter} names report both tokens of the selected variant.
    const std::string g2e5 = temp_path("resolver_snapg2e5_pw.dat");
    const std::string g2e3 = temp_path("resolver_snapg2e3_pw.dat");
    touch(g2e3);
    touch(g2e5);
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_snap_pw.dat", "_pw.dat", true,
                                      &ionic_step, &electronic_step), g2e5);
    EXPECT_EQ(ionic_step, 2);
    EXPECT_EQ(electronic_step, 5);
    std::remove(g2e3.c_str());
    std::remove(g2e5.c_str());
    // Without an electronic suffix the electronic token stays -1.
    const std::string g3 = temp_path("resolver_snapg3.txt");
    touch(g3);
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_snap.txt", ".txt", false,
                                      &ionic_step, &electronic_step), g3);
    EXPECT_EQ(ionic_step, 3);
    EXPECT_EQ(electronic_step, -1);
    std::remove(g3.c_str());
    // Null outputs keep the previous call signature working.
    EXPECT_THROW(pw::resolve_source_file(directories, "resolver_snap.txt", ".txt", false),
                 std::runtime_error);
}

TEST(PWHybridSource, ResolverFallsBackToTheSecondDirectory)
{
    const std::vector<std::string> directories = {temp_path("resolver_absent/"), testing::TempDir()};
    const std::string g1 = temp_path("resolver_psig1.txt");
    touch(g1);
    EXPECT_EQ(pw::resolve_source_file(directories, "resolver_psi.txt", ".txt", false), g1);
    std::remove(g1.c_str());
}

TEST(PWHybridSource, ResolverRejectsMissingAndElectronicSuffixWhenDisallowed)
{
    const std::vector<std::string> directories(1, testing::TempDir());
    EXPECT_THROW(pw::resolve_source_file(directories, "resolver_missing.txt", ".txt", false),
                 std::runtime_error);
    const std::string electronic = temp_path("resolver_occeig_occg1e2.txt");
    touch(electronic);
    EXPECT_THROW(pw::resolve_source_file(directories, "resolver_occeig_occ.txt", ".txt", false),
                 std::runtime_error);
    std::remove(electronic.c_str());
}
