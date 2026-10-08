#include "pw_hybrid_source.h"

#include "source_base/matrix.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <istream>
#include <stdexcept>

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

ModuleESolver::ExxSourceHeader ModuleESolver::read_exx_source_header(std::istream& in)
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

void ModuleESolver::read_exx_source_occupations(std::istream& in,
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
