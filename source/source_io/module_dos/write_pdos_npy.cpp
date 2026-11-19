#include "write_pdos_npy.h"
#include "pdos_orbitals.h"
#include "source_base/matrix.h"
#include "source_cell/unitcell.h"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <locale>
#include <random>
#include <sstream>
#include <stdexcept>
#ifdef _WIN32
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace
{
const char npy_magic[] = "\x93NUMPY\x01\x00";
const char descriptor[] = "{'descr': [('record_bytes', '<u8'),";

void write_frame(std::ostream& out, const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                 int nspin, double emin, double de, int step)
{
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "NPY PDOS requires binary64");
    const int ne = pdos.front().nc;
    std::vector<std::vector<ModuleIO::PdosOrbital>> by_type(cell.ntype);
    size_t width = 1;
    for (int it = 0; it < cell.ntype; ++it)
    {
        by_type[it] = ModuleIO::pdos_orbitals(cell.atoms[it]);
        width = std::max(width, cell.atoms[it].label.size());
    }
    width = (width + 7) / 8 * 8; // keep the float64 spectrum aligned in every row
    size_t norb = 0;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        norb += by_type[cell.iat2it[iat]].size();
    }

    // One row per (atom, l, m). The array can exceed 2 GiB without making
    // the structured dtype itself as large as the whole frame.
    const uint16_t one = 1;
    const char endian = *reinterpret_cast<const unsigned char*>(&one) == 1 ? '<' : '>';
    std::ostringstream dtype;
    dtype.imbue(std::locale::classic());
    dtype << descriptor << "('schema_version','" << endian << "i4'),('ionic_step','" << endian
          << "i4'),('nspin','" << endian << "i4'),('atom_index','" << endian << "i4'),('l','" << endian
          << "i4'),('m','" << endian << "i4'),('energy_min_eV','" << endian
          << "f8'),('energy_step_eV','" << endian << "f8'),('species','|S" << width
          << "'),('pdos','" << endian << "f8',(" << pdos.size() << "," << ne
          << "))], 'fortran_order': False, 'shape': (" << norb << ",), }";
    std::string header = dtype.str();
    header.append((64 - (10 + header.size() + 1) % 64) % 64, ' ');
    header.push_back('\n');
    const unsigned char length[] = {static_cast<unsigned char>(header.size() & 255),
                                    static_cast<unsigned char>(header.size() >> 8)};
    out.write(npy_magic, 8);
    out.write(reinterpret_cast<const char*>(length), sizeof(length));
    out.write(header.data(), header.size());

    const uint64_t row_bytes = 8 + 6 * 4 + 2 * 8 + width + 8 * static_cast<uint64_t>(ne) * pdos.size();
    if (row_bytes > std::numeric_limits<int32_t>::max())
    {
        throw std::runtime_error("one orbital exceeds the NumPy structured dtype size limit");
    }
    const uint64_t bytes = norb * row_bytes;
    unsigned char frame_size[8];
    for (unsigned int i = 0; i < 8; ++i)
    {
        frame_size[i] = static_cast<unsigned char>((bytes >> (8 * i)) & 255);
    }
    const double grid[] = {emin, de};
    std::vector<double> spectrum(ne);
    const bool spinor = nspin == 4;
    for (int iat = 0; iat < cell.nat; ++iat)
    {
        const int it = cell.iat2it[iat];
        std::string label = cell.atoms[it].label;
        label.resize(width, '\0');
        const int first = cell.get_iat2iwt()[iat];
        for (const auto& orbital : by_type[it])
        {
            const int32_t attributes[] = {1, step, nspin, iat + 1, orbital.l, orbital.m};
            out.write(reinterpret_cast<const char*>(frame_size), sizeof(frame_size));
            out.write(reinterpret_cast<const char*>(attributes), sizeof(attributes));
            out.write(reinterpret_cast<const char*>(grid), sizeof(grid));
            out.write(label.data(), width);
            for (const auto& spin : pdos)
            {
                for (int n = 0; n < ne; ++n)
                {
                    spectrum[n] = orbital.value(spin, first, n, spinor);
                }
                out.write(reinterpret_cast<const char*>(spectrum.data()), spectrum.size() * sizeof(double));
            }
        }
    }
}

void truncate_file(const std::string& path, std::streamoff size)
{
#ifdef _WIN32
    std::FILE* file = std::fopen(path.c_str(), "rb+");
    if (!file) { throw std::runtime_error(std::strerror(errno)); }
    const int status = _chsize_s(_fileno(file), size);
    std::fclose(file);
#else
    const int status = truncate(path.c_str(), size);
#endif
    if (status != 0) { throw std::runtime_error("cannot truncate an incomplete NPY record"); }
}

// Seek over the spectra. A short final record from an interrupted append can be
// discarded without rewriting the complete records, even if each energy grid differs.
std::streamoff prepare_append(const std::string& path)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in)
    {
        if (errno == ENOENT) { return 0; }
        throw std::runtime_error(std::strerror(errno));
    }
    const std::streamoff end = in.tellg();
    std::streamoff pos = 0;
    while (pos < end)
    {
        in.seekg(pos);
        unsigned char prefix[10] = {};
        in.read(reinterpret_cast<char*>(prefix), sizeof(prefix));
        const size_t got = in.gcount();
        if (in.bad() || std::memcmp(prefix, npy_magic, std::min(got, size_t(8))) != 0)
        {
            throw std::runtime_error("invalid NPY record prefix");
        }
        if (got < sizeof(prefix)) { break; }
        const int length = prefix[8] + 256 * prefix[9];
        if (length < static_cast<int>(sizeof(descriptor))) { throw std::runtime_error("invalid NPY header length"); }
        std::string header(length, '\0');
        in.read(&header[0], length);
        const size_t prefix_size = std::min<size_t>(in.gcount(), sizeof(descriptor) - 1);
        if (in.bad() || std::memcmp(header.data(), descriptor, prefix_size) != 0)
        {
            throw std::runtime_error("incompatible NPY PDOS descriptor");
        }
        if (!in) { break; }
        unsigned char raw[8];
        in.read(reinterpret_cast<char*>(raw), sizeof(raw));
        if (in.bad()) { throw std::runtime_error("error reading NPY record"); }
        if (!in) { break; }
        uint64_t bytes = 0;
        for (unsigned int i = 0; i < 8; ++i) { bytes |= uint64_t(raw[i]) << (8 * i); }
        if (bytes < 20) { throw std::runtime_error("invalid NPY payload length"); }
        const std::streamoff data_start = pos + 10 + length;
        if (bytes > static_cast<uint64_t>(end - data_start))
        {
            break;
        }
        pos = data_start + static_cast<std::streamoff>(bytes);
    }
    in.close();
    if (pos != end) { truncate_file(path, pos); }
    return pos;
}
}

std::string ModuleIO::write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                                     int nspin, double emin, double de, const std::string& out_dir,
                                     int istep, int out_app_flag)
{
    const int step = istep < 0 ? 1 : istep + 1;
    const bool numbered = out_app_flag == 2 && istep >= 0;
    const std::string suffix = numbered ? "g" + std::to_string(step) : "";
    const std::string path = out_dir + "pdos" + suffix + ".npy";
    try
    {
        const bool append = out_app_flag == 1 && istep > 0;
        const std::streamoff previous = append ? prepare_append(path) : 0;
        const std::string target = append ? path : path + ".tmp." + std::to_string(std::random_device{}());
        const auto mode = append ? std::ios::app : std::ios::trunc;
        std::ofstream out(target, std::ios::binary | mode);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try
        {
            write_frame(out, cell, pdos, nspin, emin, de, step);
            out.close();
            if (!append)
            {
#ifdef _WIN32
                if (!MoveFileExA(target.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING))
                {
                    throw std::runtime_error("cannot replace NPY file (Windows error "
                                             + std::to_string(GetLastError()) + ")");
                }
#else
                if (std::rename(target.c_str(), path.c_str()) != 0)
                {
                    throw std::runtime_error(std::strerror(errno));
                }
#endif
            }
        }
        catch (...)
        {
            out.exceptions(std::ios::goodbit);
            if (out.is_open()) { out.close(); }
            if (append) { truncate_file(path, previous); }
            else { std::remove(target.c_str()); }
            throw;
        }
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(path + ": " + e.what());
    }
    return path;
}
