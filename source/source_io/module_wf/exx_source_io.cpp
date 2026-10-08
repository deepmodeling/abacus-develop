#include "exx_source_io.h"

#include "source_base/matrix.h"
#include "source_cell/klist.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <istream>
#include <ostream>
#include <stdexcept>

void ModuleIO::write_exx_source(std::ostream& out,
                              const K_Vectors& points,
                              const ModuleBase::matrix& weights,
                              int nspin,
                              int nx,
                              int ny,
                              int nz,
                              const std::string& configuration)
{
    out << "ABACUS_EXX_SOURCE_V1\n" << configuration << '\n';
    out << nspin << ' ' << nx << ' ' << ny << ' ' << nz << ' ' << weights.nr << ' ' << weights.nc << '\n';
    out << points.nmp[0] << ' ' << points.nmp[1] << ' ' << points.nmp[2] << '\n';
    out << std::setprecision(17);
    for (int iq = 0; iq < weights.nr; ++iq)
    {
        const auto& q = points.kvec_d[iq];
        out << q.x << ' ' << q.y << ' ' << q.z << ' ' << points.wk[iq] << ' ' << points.isk[iq];
        for (int ib = 0; ib < weights.nc; ++ib)
        {
            out << ' ' << weights(iq, ib);
        }
        out << '\n';
    }
    if (!out)
    {
        throw std::runtime_error("Cannot write EXX source checkpoint");
    }
}

void ModuleIO::read_exx_source(std::istream& in,
                             K_Vectors& points,
                             ModuleBase::matrix& weights,
                             int nspin,
                             int nx,
                             int ny,
                             int nz,
                             const std::string& configuration)
{
    std::string magic;
    std::string saved_configuration;
    std::getline(in, magic);
    std::getline(in, saved_configuration);
    int saved_spin = 0;
    int saved_nx = 0;
    int saved_ny = 0;
    int saved_nz = 0;
    int nks = 0;
    int nbands = 0;
    in >> saved_spin >> saved_nx >> saved_ny >> saved_nz >> nks >> nbands;
    if (!in || magic != "ABACUS_EXX_SOURCE_V1" || saved_configuration != configuration
        || saved_spin != nspin || saved_nx != nx || saved_ny != ny || saved_nz != nz
        || nks <= 0 || nbands <= 0 || nks > 100000 || nbands > 100000
        || static_cast<long long>(nks) * nbands > 100000000)
    {
        throw std::runtime_error("Missing, incompatible or invalid EXX source checkpoint; regenerate with a converged hybrid SCF");
    }
    const int spin_mult = nspin == 2 ? 2 : 1;
    if (nks % spin_mult != 0)
    {
        throw std::runtime_error("Invalid EXX source spin dimensions");
    }
    in >> points.nmp[0] >> points.nmp[1] >> points.nmp[2];
    points.set_nks(nks);
    points.set_nkstot(nks);
    const int nq = nks / spin_mult;
    points.set_nkstot_nospin(nq);
    points.set_spin_mult(spin_mult);
    points.kvec_d.resize(nks);
    points.kvec_c.resize(nks);
    points.wk.resize(nks);
    points.isk.resize(nks);
    points.ngk.resize(nks);
    points.ik2iktot.resize(nks);
    weights.create(nks, nbands);
    const double expected_weight = nspin == 1 ? 2.0 / nq : 1.0 / nq;
    for (int iq = 0; iq < nks; ++iq)
    {
        auto& q = points.kvec_d[iq];
        in >> q.x >> q.y >> q.z >> points.wk[iq] >> points.isk[iq];
        if (!in || !std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z)
            || !std::isfinite(points.wk[iq]) || std::abs(points.wk[iq] - expected_weight) > 1e-10
            || points.isk[iq] != iq / nq)
        {
            throw std::runtime_error("EXX source requires a complete uniform q mesh with ordered spin channels");
        }
        points.ik2iktot[iq] = iq;
        for (int ib = 0; ib < nbands; ++ib)
        {
            in >> weights(iq, ib);
            if (!in || !std::isfinite(weights(iq, ib)) || weights(iq, ib) < 0
                || weights(iq, ib) > expected_weight + 1e-10)
            {
                throw std::runtime_error("Invalid EXX source occupation");
            }
        }
    }
    in >> std::ws;
    if (!in.eof())
    {
        throw std::runtime_error("Unexpected trailing EXX source checkpoint data");
    }
}

namespace
{
template <typename T>
void read_binary(std::istream& in, T& value)
{
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in)
    {
        throw std::runtime_error("Missing or truncated EXX source wavefunction");
    }
}
}

ModuleIO::ExxSourceHeader ModuleIO::read_exx_source_header(std::istream& in)
{
    ExxSourceHeader header{};
    int record = 0;
    read_binary(in, record);
    if (record != 72)
    {
        throw std::runtime_error("Invalid EXX source wavefunction header");
    }
    read_binary(in, header.ik);
    read_binary(in, header.nks);
    for (double& coordinate : header.kvec_c)
    {
        read_binary(in, coordinate);
        if (!std::isfinite(coordinate))
        {
            throw std::runtime_error("Invalid EXX source wavefunction k point");
        }
    }
    read_binary(in, header.weight);
    read_binary(in, header.npw);
    read_binary(in, header.nbands);
    read_binary(in, header.ecutwfc);
    read_binary(in, header.lat0);
    read_binary(in, header.tpiba);
    read_binary(in, record);
    if (record != 72 || header.ik <= 0 || header.nks < header.ik || header.nks > 100000
        || header.nbands <= 0 || header.nbands > 100000 || header.npw <= 0
        || !std::isfinite(header.weight) || header.weight <= 0
        || !std::isfinite(header.ecutwfc) || header.ecutwfc <= 0
        || !std::isfinite(header.lat0) || header.lat0 <= 0
        || !std::isfinite(header.tpiba) || header.tpiba <= 0)
    {
        throw std::runtime_error("Invalid EXX source wavefunction dimensions or metadata");
    }
    // Two 72-byte records, a Miller-index record, and complex<double> band records.
    const long long expected_bytes = 160LL + 8 + 12LL * header.npw
                                    + header.nbands * (8LL + 16LL * header.npw);
    in.seekg(0, std::ios::end);
    if (!in || in.tellg() != expected_bytes)
    {
        throw std::runtime_error("Truncated or invalid EXX source wavefunction size");
    }
    return header;
}

void ModuleIO::read_exx_source_occupations(std::istream& in,
                                          const std::vector<ExxSourceHeader>& headers,
                                          int nspin,
                                          ModuleBase::matrix& weights)
{
    const int spin_mult = nspin == 2 ? 2 : 1;
    const int nks = headers.size();
    if ((nspin != 1 && nspin != 2) || nks == 0 || nks % spin_mult != 0)
    {
        throw std::runtime_error("Invalid EXX source spin dimensions");
    }
    const int nq = nks / spin_mult;
    const int nbands = headers.front().nbands;
    if (static_cast<long long>(nks) * nbands > 100000000)
    {
        throw std::runtime_error("Invalid EXX source occupation dimensions");
    }
    std::string line;
    std::getline(in, line);
    std::getline(in, line);
    std::getline(in, line);
    int saved_spin = 0;
    if (std::sscanf(line.c_str(), " Spin number %d", &saved_spin) != 1 || saved_spin != nspin)
    {
        throw std::runtime_error("Missing or incompatible EXX source eig_occ.txt");
    }
    weights.create(nks, nbands);
    const double expected_weight = nspin == 1 ? 2.0 / nq : 1.0 / nq;
    for (int iq = 0; iq < nks; ++iq)
    {
        const auto& header = headers[iq];
        do
        {
            if (!std::getline(in, line))
            {
                throw std::runtime_error("Truncated EXX source eig_occ.txt");
            }
        } while (line.find_first_not_of(" \t\r") == std::string::npos);
        int spin = 0;
        int index = 0;
        int total = 0;
        int npw = 0;
        double q[3] = {};
        const int fields = std::sscanf(line.c_str(),
            " spin=%d k-point=%d/%d Cartesian=%lf %lf %lf (%d plane wave)",
            &spin, &index, &total, &q[0], &q[1], &q[2], &npw);
        if (fields != 7 || spin != iq / nq + 1 || index != iq % nq + 1 || total != nq
            || npw != header.npw || header.ik != iq + 1 || header.nks != nks
            || header.nbands != nbands || std::abs(header.weight - expected_weight) > 1e-10)
        {
            throw std::runtime_error("EXX source requires matching eig_occ.txt and a complete uniform q mesh");
        }
        for (int axis = 0; axis < 3; ++axis)
        {
            // eig_occ.txt prints Cartesian coordinates with at least eight significant digits.
            const double tolerance = 5e-8 * std::max(1.0, std::abs(header.kvec_c[axis]));
            if (!std::isfinite(q[axis]) || std::abs(q[axis] - header.kvec_c[axis]) > tolerance)
            {
                throw std::runtime_error("EXX source eig_occ.txt k point differs from wavefunction");
            }
        }
        for (int ib = 0; ib < nbands; ++ib)
        {
            std::getline(in, line);
            std::istringstream band(line);
            double energy = 0;
            band >> index >> energy >> weights(iq, ib);
            const bool parsed = static_cast<bool>(band);
            std::string extra;
            const bool trailing = static_cast<bool>(band >> extra);
            if (!parsed || trailing || index != ib + 1 || !std::isfinite(energy)
                || !std::isfinite(weights(iq, ib)) || weights(iq, ib) < 0
                || weights(iq, ib) > expected_weight + 1e-10)
            {
                throw std::runtime_error("Invalid EXX source occupation in eig_occ.txt");
            }
        }
    }
    in >> std::ws;
    if (!in.eof())
    {
        throw std::runtime_error("Unexpected trailing EXX source eig_occ.txt data; use a single SCF snapshot");
    }
}
