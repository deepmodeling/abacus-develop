#include "pw_hybrid_source.h"

#include "source_base/matrix.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
template <typename T>
void read_binary(std::istream &in, T &value)
{
    in.read(reinterpret_cast<char *>(&value), sizeof(T));
    if (!in)
    {
        throw std::runtime_error("Missing or truncated EXX source wavefunction");
    }
}
} // namespace

pw::ExxSourceHeader pw::read_exx_source_header(std::istream &in)
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
    for (double &coordinate : header.kvec_c)
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
    if (record != 72 || header.ik <= 0 || header.nks < header.ik || header.nks > 100000 || header.nbands <= 0 ||
        header.nbands > 100000 || header.npw <= 0 || !std::isfinite(header.weight) || header.weight <= 0 ||
        !std::isfinite(header.ecutwfc) || header.ecutwfc <= 0 || !std::isfinite(header.lat0) || header.lat0 <= 0 ||
        !std::isfinite(header.tpiba) || header.tpiba <= 0)
    {
        throw std::runtime_error("Invalid EXX source wavefunction dimensions or metadata");
    }
    // Two 72-byte records, a Miller-index record, and complex<double> band records.
    const long long expected_bytes = 160LL + 8 + 12LL * header.npw + header.nbands * (8LL + 16LL * header.npw);
    in.seekg(0, std::ios::end);
    if (!in || in.tellg() != expected_bytes)
    {
        throw std::runtime_error("Truncated or invalid EXX source wavefunction size");
    }
    return header;
}

void pw::read_exx_source_occupations(
    std::istream &in,
    const std::vector<ExxSourceHeader> &headers,
    int nspin,
    ModuleBase::matrix &weights)
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
    int saved_spin = 0;
    bool spin_found = false;
    // The header may carry provenance comments and the ionic-step marker before
    // the "Spin number" line, so scan until the spin line is located.
    while (std::getline(in, line))
    {
        if (std::sscanf(line.c_str(), " Spin number %d", &saved_spin) == 1)
        {
            spin_found = true;
            break;
        }
        if (line.find(" spin=") != std::string::npos)
        {
            break;
        }
    }
    if (!spin_found || saved_spin != nspin)
    {
        throw std::runtime_error("Missing or incompatible EXX source eig_occ.txt");
    }
    weights.create(nks, nbands);
    const double expected_weight = nspin == 1 ? 2.0 / nq : 1.0 / nq;
    for (int iq = 0; iq < nks; ++iq)
    {
        const pw::ExxSourceHeader &header = headers[iq];
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
        const int fields = std::sscanf(
            line.c_str(), " spin=%d k-point=%d/%d Cartesian=%lf %lf %lf (%d plane wave)", &spin, &index, &total, &q[0],
            &q[1], &q[2], &npw);
        if (fields != 7 || spin != iq / nq + 1 || index != iq % nq + 1 || total != nq || npw != header.npw ||
            header.ik != iq + 1 || header.nks != nks || header.nbands != nbands ||
            std::abs(header.weight - expected_weight) > 1e-10)
        {
            throw std::runtime_error("EXX source requires matching eig_occ.txt and a complete uniform q mesh");
        }
        for (int axis = 0; axis < 3; ++axis)
        {
            if (nspin == 2 && iq >= nq && header.kvec_c[axis] != headers[iq - nq].kvec_c[axis])
            {
                throw std::runtime_error("EXX source spin channels must use the same q mesh");
            }
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
            if (!parsed || trailing || index != ib + 1 || !std::isfinite(energy) || !std::isfinite(weights(iq, ib)))
            {
                throw std::runtime_error("Invalid EXX source occupation in eig_occ.txt");
            }
            if (weights(iq, ib) < 0 || weights(iq, ib) > expected_weight + 1e-10)
            {
                // Methfessel-Paxton smearing can write occupations below zero or above
                // the k-point weight, and the EXX operator silently drops wg < 1e-12;
                // only non-negative occupation schemes are supported.
                throw std::runtime_error(
                    "EXX source eig_occ.txt occupation is outside [0, k-point weight]; "
                    "Methfessel-Paxton smearing (mp/mp2/mp3) can produce such occupations, "
                    "which hybrid NSCF does not support. Rerun the SCF source with "
                    "smearing_method fixed, gauss, mv or fd");
            }
        }
    }
    in >> std::ws;
    if (!in.eof())
    {
        throw std::runtime_error("Unexpected trailing EXX source eig_occ.txt data; use a single SCF snapshot");
    }
}

namespace
{
/**
 * @brief Parse the optional "g{step}" and "e{iter}" tokens between the file
 *        stem and the fixed suffix.
 * @param tokens Candidate token string, e.g. "g2" or "g2e17".
 * @param ionic_step Parsed one-based ionic step, or -1.
 * @param electronic_step Parsed one-based electronic iteration, or -1.
 * @return True when the token string is a valid g/e combination.
 */
bool parse_step_tokens(const std::string &tokens, int &ionic_step, int &electronic_step)
{
    ionic_step = -1;
    electronic_step = -1;
    if (tokens.empty() || tokens[0] != 'g')
    {
        return false;
    }
    size_t position = 1;
    const size_t ionic_begin = position;
    while (position < tokens.size() && std::isdigit(static_cast<unsigned char>(tokens[position])) != 0)
    {
        ++position;
    }
    if (position == ionic_begin)
    {
        return false;
    }
    ionic_step = std::stoi(tokens.substr(ionic_begin, position - ionic_begin));
    if (position == tokens.size())
    {
        return true;
    }
    if (tokens[position] != 'e')
    {
        return false;
    }
    ++position;
    const size_t electronic_begin = position;
    while (position < tokens.size() && std::isdigit(static_cast<unsigned char>(tokens[position])) != 0)
    {
        ++position;
    }
    if (position == electronic_begin || position != tokens.size())
    {
        return false;
    }
    electronic_step = std::stoi(tokens.substr(electronic_begin, position - electronic_begin));
    return true;
}

/**
 * @brief Find one source file inside a single directory.
 * @return Full path, or an empty string when no candidate exists.
 */
std::string find_in_directory(
    const std::string &directory,
    const std::string &fixed_name,
    const std::string &tail,
    const bool allow_electronic_suffix,
    int *ionic_step,
    int *electronic_step)
{
    const std::string fixed_path = directory + fixed_name;
    std::ifstream fixed_stream(fixed_path.c_str(), std::ios::binary);
    if (fixed_stream.is_open())
    {
        if (ionic_step != nullptr)
        {
            *ionic_step = -1;
        }
        if (electronic_step != nullptr)
        {
            *electronic_step = -1;
        }
        return fixed_path;
    }

    const std::string stem = fixed_name.substr(0, fixed_name.size() - tail.size());
    DIR *directory_stream = opendir(directory.c_str());
    if (directory_stream == nullptr)
    {
        return std::string();
    }

    std::string best_name;
    int best_ionic_step = -1;
    int best_electronic_step = -1;
    while (const dirent *entry = readdir(directory_stream))
    {
        const std::string name(entry->d_name);
        const size_t minimum_size = stem.size() + tail.size() + 2; // at least "g1"
        if (name.size() < minimum_size)
        {
            continue;
        }
        if (name.compare(0, stem.size(), stem) != 0)
        {
            continue;
        }
        if (name.compare(name.size() - tail.size(), tail.size(), tail) != 0)
        {
            continue;
        }
        const std::string tokens = name.substr(stem.size(), name.size() - stem.size() - tail.size());
        int ionic_step = -1;
        int electronic_step = -1;
        if (!parse_step_tokens(tokens, ionic_step, electronic_step))
        {
            continue;
        }
        if (ionic_step < 1)
        {
            continue;
        }
        if (!allow_electronic_suffix && electronic_step >= 0)
        {
            continue;
        }
        if (ionic_step > best_ionic_step || (ionic_step == best_ionic_step && electronic_step > best_electronic_step))
        {
            best_ionic_step = ionic_step;
            best_electronic_step = electronic_step;
            best_name = name;
        }
    }
    closedir(directory_stream);

    if (best_name.empty())
    {
        return std::string();
    }
    if (ionic_step != nullptr)
    {
        *ionic_step = best_ionic_step;
    }
    if (electronic_step != nullptr)
    {
        *electronic_step = best_electronic_step;
    }
    return directory + best_name;
}
} // namespace

std::string pw::resolve_source_file(
    const std::vector<std::string> &directories,
    const std::string &fixed_name,
    const std::string &tail,
    const bool allow_electronic_suffix,
    int *ionic_step,
    int *electronic_step)
{
    if (directories.empty() || fixed_name.empty() || tail.empty())
    {
        throw std::runtime_error("Invalid EXX source file request");
    }
    if (fixed_name.size() <= tail.size() || fixed_name.compare(fixed_name.size() - tail.size(), tail.size(), tail) != 0)
    {
        throw std::runtime_error("EXX source fixed name does not end with the expected suffix");
    }

    for (const std::string &directory : directories)
    {
        if (directory.empty())
        {
            continue;
        }
        const std::string found =
            find_in_directory(directory, fixed_name, tail, allow_electronic_suffix, ionic_step, electronic_step);
        if (!found.empty())
        {
            return found;
        }
    }

    throw std::runtime_error(
        "Missing or incompatible EXX source file: " + fixed_name + " (expected the fixed name or a g{step} variant)");
}
