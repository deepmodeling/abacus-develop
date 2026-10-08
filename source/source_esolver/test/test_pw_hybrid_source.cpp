#include "source_esolver/pw_hybrid_source.h"
#include "source_base/matrix.h"

#include "gtest/gtest.h"

#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace
{
std::vector<ModuleESolver::ExxSourceHeader> source_headers(int nspin)
{
    const int nks = nspin == 2 ? 4 : 2;
    std::vector<ModuleESolver::ExxSourceHeader> headers(nks);
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
    const auto headers = source_headers(nspin);
    std::ostringstream out;
    out << "1 # ionic step\n Electronic state energy (eV) and occupations\n Spin number " << nspin << '\n';
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
        const auto headers = source_headers(nspin);
        std::istringstream in(eig_occ(nspin));
        ModuleBase::matrix weights;
        ModuleESolver::read_exx_source_occupations(in, headers, nspin, weights);
        EXPECT_EQ(weights.nr, headers.size());
        EXPECT_EQ(weights.nc, 3);
        EXPECT_DOUBLE_EQ(weights(1, 1), headers[1].weight * 0.125);
    }
}

TEST(PWHybridSource, RejectsInvalidExistingOccupationFiles)
{
    const auto headers = source_headers(1);
    ModuleBase::matrix weights;
    const std::string saved = eig_occ(1);
    std::istringstream truncated(saved.substr(0, saved.size() - 6));
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(truncated, headers, 1, weights), std::runtime_error);
    std::istringstream trailing(saved + saved);
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(trailing, headers, 1, weights), std::runtime_error);
    std::istringstream wrong_spin(saved);
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(wrong_spin, headers, 2, weights), std::runtime_error);
    auto changed = headers;
    changed[1].kvec_c[0] += 0.1;
    std::istringstream wrong_coordinates(saved);
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(wrong_coordinates, changed, 1, weights), std::runtime_error);
    std::string invalid = saved;
    invalid.replace(invalid.find("0.125"), 5, "-0.25");
    std::istringstream negative(invalid);
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(negative, headers, 1, weights), std::runtime_error);
}

TEST(PWHybridSource, RejectsDifferentMeshesInSpinChannels)
{
    auto headers = source_headers(2);
    headers[2].kvec_c[0] = 0.2;
    std::string saved = eig_occ(2);
    const std::string original = "spin=2 k-point=1/2 Cartesian=0";
    const std::size_t location = saved.find(original);
    ASSERT_NE(location, std::string::npos);
    saved.replace(location, original.size(), "spin=2 k-point=1/2 Cartesian=0.2");
    std::istringstream in(saved);
    ModuleBase::matrix weights;
    EXPECT_THROW(ModuleESolver::read_exx_source_occupations(in, headers, 2, weights), std::runtime_error);
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
    const auto header = ModuleESolver::read_exx_source_header(in);
    EXPECT_EQ(header.nbands, 3);
    EXPECT_DOUBLE_EQ(header.kvec_c[0], coordinate);
    EXPECT_DOUBLE_EQ(header.weight, weight);
    std::istringstream truncated(saved.substr(0, saved.size() - 1), std::ios::binary);
    EXPECT_THROW(ModuleESolver::read_exx_source_header(truncated), std::runtime_error);
    saved[0] = 0;
    std::istringstream invalid(saved, std::ios::binary);
    EXPECT_THROW(ModuleESolver::read_exx_source_header(invalid), std::runtime_error);
}
