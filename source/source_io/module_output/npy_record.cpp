#include "npy_record.h"

#include <algorithm>
#include <cerrno>
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
std::streamoff prepare_append(const std::string& path, const std::string& identifying_field)
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
            throw std::runtime_error("incompatible NPY record descriptor");
        }
        if (!in) { break; }
        if (header.find("('" + identifying_field + "',") == std::string::npos)
        {
            throw std::runtime_error("unexpected dataset in NPY file");
        }
        unsigned char raw[8];
        in.read(reinterpret_cast<char*>(raw), sizeof(raw));
        if (in.bad()) { throw std::runtime_error("error reading NPY record"); }
        if (!in) { break; }
        uint64_t bytes = 0;
        for (unsigned int i = 0; i < 8; ++i) { bytes |= uint64_t(raw[i]) << (8 * i); }
        if (bytes < sizeof(ModuleIO::NpyPrefix)) { throw std::runtime_error("invalid NPY payload length"); }
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

void ModuleIO::write_npy_record(const std::string& path, bool append, int step,
                                const std::string& fields, const std::string& identifying_field,
                                std::uint64_t count, std::uint64_t row_bytes, const NpyRows& write_rows)
{
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "NPY requires binary64");
    if (count == 0 || row_bytes < 16 || row_bytes > std::numeric_limits<int32_t>::max()
        || count > static_cast<uint64_t>(std::numeric_limits<std::streamoff>::max()) / row_bytes)
    {
        throw std::runtime_error(path + ": NPY record dimensions exceed the structured dtype or file size limit");
    }
    const std::uint16_t one = 1;
    const char endian = *reinterpret_cast<const unsigned char*>(&one) == 1 ? '<' : '>';
    std::string native_fields = fields;
    std::replace(native_fields.begin(), native_fields.end(), '=', endian);
    std::ostringstream dtype;
    dtype.imbue(std::locale::classic());
    dtype << descriptor << "('schema_version','" << endian << "i4'),('ionic_step','" << endian << "i4'),"
          << native_fields << "], 'fortran_order': False, 'shape': (" << count << ",), }";
    std::string header = dtype.str();
    header.append((64 - (10 + header.size() + 1) % 64) % 64, ' ');
    header.push_back('\n');
    if (header.size() > 65535) { throw std::runtime_error(path + ": NPY v1 header is too long"); }
    const unsigned char length[] = {static_cast<unsigned char>(header.size() & 255),
                                    static_cast<unsigned char>(header.size() >> 8)};
    NpyPrefix prefix;
    const uint64_t bytes = count * row_bytes;
    for (unsigned int i = 0; i < 8; ++i) { prefix[i] = static_cast<char>((bytes >> (8 * i)) & 255); }
    const std::int32_t attributes[] = {1, step};
    std::memcpy(prefix.data() + 8, attributes, sizeof(attributes));
    try
    {
        const std::streamoff previous = append ? prepare_append(path, identifying_field) : 0;
        const std::string target = append ? path : path + ".tmp." + std::to_string(std::random_device{}());
        const auto mode = append ? std::ios::app : std::ios::trunc;
        std::ofstream out(target, std::ios::binary | mode);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try
        {
            out.write(npy_magic, 8);
            out.write(reinterpret_cast<const char*>(length), sizeof(length));
            out.write(header.data(), header.size());
            const auto start = out.tellp();
            write_rows(out, prefix);
            if (out.tellp() - start != static_cast<std::streamoff>(bytes))
            {
                throw std::runtime_error("NPY payload does not match the declared dimensions");
            }
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
}
