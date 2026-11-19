#include "npy_record.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace
{
const char magic[] = "\x93NUMPY\x01\x00";
const std::string descriptor = "{'descr': [('record_bytes', '<u8'),";

// Only the header and first eight payload bytes of each previous step are read.
void check_append(const std::string& path, const char* field)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in && errno == ENOENT) { return; }
    in.exceptions(std::ios::badbit | std::ios::failbit);
    const std::streamoff end = in.tellg();
    in.seekg(0);
    while (in.tellg() < end)
    {
        unsigned char tag[10];
        in.read(reinterpret_cast<char*>(tag), sizeof(tag));
        const int length = tag[8] + 256 * tag[9];
        if (std::memcmp(tag, magic, 8) != 0 || length < static_cast<int>(descriptor.size()))
        {
            throw std::runtime_error("invalid NPY header");
        }
        std::string header(length, '\0');
        in.read(&header[0], length);
        if (header.compare(0, descriptor.size(), descriptor) != 0
            || header.find(std::string("('") + field + "',") == std::string::npos)
        {
            throw std::runtime_error("incompatible DOS array");
        }
        const std::streamoff start = in.tellg();
        unsigned char size[8];
        in.read(reinterpret_cast<char*>(size), sizeof(size));
        std::uint64_t bytes = 0;
        for (int i = 7; i >= 0; --i) { bytes = (bytes << 8) | size[i]; }
        if (bytes < 16 || bytes > static_cast<std::uint64_t>(end - start))
        {
            throw std::runtime_error("incomplete DOS array");
        }
        in.seekg(start + static_cast<std::streamoff>(bytes));
    }
}
}

std::array<char, 16> ModuleIO::open_npy_record(std::ofstream& out, const std::string& path, bool append, int step,
                                              const std::string& fields, const char* field, std::uint64_t count,
                                              std::uint64_t row_bytes)
{
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "NPY requires binary64");
    if (count == 0 || row_bytes < 16 || row_bytes > std::numeric_limits<std::int32_t>::max()
        || count > static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()) / row_bytes)
    {
        throw std::runtime_error(path + ": unsupported DOS array dimensions");
    }
    const std::uint16_t one = 1;
    const char endian = *reinterpret_cast<const unsigned char*>(&one) == 1 ? '<' : '>';
    std::string header = descriptor + "('schema_version','=i4'),('ionic_step','=i4')," + fields
                       + "], 'fortran_order': False, 'shape': (" + std::to_string(count) + ",), }";
    std::replace(header.begin(), header.end(), '=', endian);
    header.append((64 - (11 + header.size()) % 64) % 64, ' ');
    header += '\n';
    if (header.size() > 65535) { throw std::runtime_error(path + ": NPY v1 header is too long"); }
    const unsigned char length[] = {static_cast<unsigned char>(header.size()), static_cast<unsigned char>(header.size() >> 8)};
    try
    {
        if (append) { check_append(path, field); }
        out.exceptions(std::ios::badbit | std::ios::failbit);
        out.open(path, std::ios::binary | (append ? std::ios::app : std::ios::trunc));
        out.write(magic, 8);
        out.write(reinterpret_cast<const char*>(length), sizeof(length));
        out.write(header.data(), header.size());
    }
    catch (const std::exception& e) { throw std::runtime_error(path + ": " + e.what()); }

    std::array<char, 16> prefix;
    const std::uint64_t bytes = count * row_bytes;
    for (int i = 0; i < 8; ++i) { prefix[i] = static_cast<char>(bytes >> (8 * i)); }
    const std::int32_t attributes[] = {1, step};
    std::memcpy(prefix.data() + 8, attributes, sizeof(attributes));
    return prefix;
}
