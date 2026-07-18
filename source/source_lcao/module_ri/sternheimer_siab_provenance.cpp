#include "sternheimer_siab_provenance.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

#if defined(__linux__)
#include <unistd.h>
#endif

namespace module_ri
{
namespace sternheimer_siab
{
namespace
{

constexpr std::array<std::uint32_t, 64> kRoundConstants{{
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
    0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
    0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
    0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
    0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
    0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
    0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
    0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U, 0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U,
}};

std::uint32_t rotate_right(const std::uint32_t value, const unsigned int shift)
{
    return (value >> shift) | (value << (32U - shift));
}

bool is_hex_string(const std::string& value, const std::size_t size)
{
    return value.size() == size && std::all_of(value.begin(), value.end(), [](const unsigned char character) {
               return (character >= '0' && character <= '9') || (character >= 'a' && character <= 'f')
                      || (character >= 'A' && character <= 'F');
           });
}

std::string normalized_commit_token(const std::string& value)
{
    std::size_t begin = 0;
    while (begin != value.size() && std::isspace(static_cast<unsigned char>(value[begin])))
    {
        ++begin;
    }
    std::size_t end = begin;
    while (end != value.size() && !std::isspace(static_cast<unsigned char>(value[end])) && value[end] != '(')
    {
        ++end;
    }
    return value.substr(begin, end - begin);
}

std::string join_path(const std::string& directory, const std::string& filename)
{
    if (filename.empty() || (!filename.empty() && filename.front() == '/') || directory.empty())
    {
        return filename;
    }
    return directory.back() == '/' ? directory + filename : directory + "/" + filename;
}

void update_u64(Sha256& digest, const std::uint64_t value)
{
    std::array<unsigned char, 8> bytes;
    for (std::size_t index = 0; index != bytes.size(); ++index)
    {
        bytes[index] = static_cast<unsigned char>((value >> (56U - 8U * index)) & 0xffU);
    }
    digest.update(bytes.data(), bytes.size());
}

void update_file(Sha256& digest, const std::string& path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("Cannot open required Sternheimer SIAB provenance file: " + path);
    }
    std::array<char, 1024 * 1024> buffer;
    while (input)
    {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize count = input.gcount();
        if (count > 0)
        {
            digest.update(reinterpret_cast<const unsigned char*>(buffer.data()), static_cast<std::size_t>(count));
        }
    }
    if (!input.eof())
    {
        throw std::runtime_error("Failed while hashing Sternheimer SIAB provenance file: " + path);
    }
}

std::uint64_t file_size(const std::string& path)
{
    std::ifstream input(path.c_str(), std::ios::binary | std::ios::ate);
    if (!input)
    {
        throw std::runtime_error("Cannot open required Sternheimer SIAB provenance file: " + path);
    }
    const std::streamoff position = input.tellg();
    if (position < 0)
    {
        throw std::runtime_error("Cannot determine size of Sternheimer SIAB provenance file: " + path);
    }
    return static_cast<std::uint64_t>(position);
}

} // namespace

Sha256::Sha256()
    : state_{{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU, 0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}}
{
    buffer_.fill(0);
}

void Sha256::transform(const unsigned char* block)
{
    std::array<std::uint32_t, 64> words;
    for (std::size_t index = 0; index != 16; ++index)
    {
        words[index] = (static_cast<std::uint32_t>(block[4 * index]) << 24U)
                       | (static_cast<std::uint32_t>(block[4 * index + 1]) << 16U)
                       | (static_cast<std::uint32_t>(block[4 * index + 2]) << 8U)
                       | static_cast<std::uint32_t>(block[4 * index + 3]);
    }
    for (std::size_t index = 16; index != words.size(); ++index)
    {
        const std::uint32_t s0
            = rotate_right(words[index - 15], 7) ^ rotate_right(words[index - 15], 18) ^ (words[index - 15] >> 3U);
        const std::uint32_t s1
            = rotate_right(words[index - 2], 17) ^ rotate_right(words[index - 2], 19) ^ (words[index - 2] >> 10U);
        words[index] = words[index - 16] + s0 + words[index - 7] + s1;
    }

    std::uint32_t a = state_[0];
    std::uint32_t b = state_[1];
    std::uint32_t c = state_[2];
    std::uint32_t d = state_[3];
    std::uint32_t e = state_[4];
    std::uint32_t f = state_[5];
    std::uint32_t g = state_[6];
    std::uint32_t h = state_[7];
    for (std::size_t index = 0; index != words.size(); ++index)
    {
        const std::uint32_t sigma1 = rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
        const std::uint32_t choice = (e & f) ^ (~e & g);
        const std::uint32_t temporary1 = h + sigma1 + choice + kRoundConstants[index] + words[index];
        const std::uint32_t sigma0 = rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
        const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
        const std::uint32_t temporary2 = sigma0 + majority;
        h = g;
        g = f;
        f = e;
        e = d + temporary1;
        d = c;
        c = b;
        b = a;
        a = temporary1 + temporary2;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
}

void Sha256::update(const unsigned char* data, const std::size_t size)
{
    if (finished_)
    {
        throw std::logic_error("Cannot update a finalized SHA256 digest.");
    }
    if (size != 0 && data == nullptr)
    {
        throw std::invalid_argument("Cannot hash a null non-empty byte range.");
    }
    if (size > std::numeric_limits<std::uint64_t>::max() - total_size_)
    {
        throw std::overflow_error("SHA256 input length exceeds uint64 capacity.");
    }
    total_size_ += static_cast<std::uint64_t>(size);
    std::size_t consumed = 0;
    if (buffer_size_ != 0)
    {
        const std::size_t take = std::min(size, buffer_.size() - buffer_size_);
        std::memcpy(buffer_.data() + buffer_size_, data, take);
        buffer_size_ += take;
        consumed += take;
        if (buffer_size_ == buffer_.size())
        {
            transform(buffer_.data());
            buffer_size_ = 0;
        }
    }
    while (size - consumed >= buffer_.size())
    {
        transform(data + consumed);
        consumed += buffer_.size();
    }
    if (consumed != size)
    {
        buffer_size_ = size - consumed;
        std::memcpy(buffer_.data(), data + consumed, buffer_size_);
    }
}

std::string Sha256::finish()
{
    if (finished_)
    {
        throw std::logic_error("Cannot finalize a SHA256 digest twice.");
    }
    if (total_size_ > std::numeric_limits<std::uint64_t>::max() / 8U)
    {
        throw std::overflow_error("SHA256 bit length exceeds uint64 capacity.");
    }
    const std::uint64_t bit_size = total_size_ * 8U;
    buffer_[buffer_size_++] = 0x80U;
    if (buffer_size_ > 56)
    {
        std::fill(buffer_.begin() + static_cast<std::ptrdiff_t>(buffer_size_), buffer_.end(), 0U);
        transform(buffer_.data());
        buffer_size_ = 0;
    }
    std::fill(buffer_.begin() + static_cast<std::ptrdiff_t>(buffer_size_), buffer_.begin() + 56, 0U);
    for (std::size_t index = 0; index != 8; ++index)
    {
        buffer_[56 + index] = static_cast<unsigned char>((bit_size >> (56U - 8U * index)) & 0xffU);
    }
    transform(buffer_.data());
    finished_ = true;

    std::ostringstream output;
    output << std::hex << std::nouppercase << std::setfill('0');
    for (const std::uint32_t word: state_)
    {
        output << std::setw(8) << word;
    }
    return output.str();
}

std::string sha256_bytes(const std::vector<unsigned char>& bytes)
{
    Sha256 digest;
    digest.update(bytes.data(), bytes.size());
    return digest.finish();
}

std::string sha256_file(const std::string& path)
{
    Sha256 digest;
    update_file(digest, path);
    return digest.finish();
}

std::string sha256_file_manifest(const std::vector<std::string>& paths)
{
    if (paths.empty())
    {
        throw std::invalid_argument("Sternheimer SIAB provenance file manifest must not be empty.");
    }
    if (paths.size() == 1)
    {
        return sha256_file(paths.front());
    }

    Sha256 digest;
    const std::string marker = "ABACUS_STERNHEIMER_SIAB_FILE_MANIFEST_V1";
    digest.update(reinterpret_cast<const unsigned char*>(marker.data()), marker.size());
    update_u64(digest, static_cast<std::uint64_t>(paths.size()));
    for (const std::string& path: paths)
    {
        update_u64(digest, file_size(path));
        update_file(digest, path);
    }
    return digest.finish();
}

std::vector<std::string> resolve_required_input_files(const std::string& directory,
                                                      const std::vector<std::string>& filenames,
                                                      const std::string& label)
{
    if (filenames.empty())
    {
        throw std::runtime_error("Sternheimer SIAB production provenance requires at least one " + label + " file.");
    }
    std::vector<std::string> paths;
    paths.reserve(filenames.size());
    for (const std::string& filename: filenames)
    {
        if (filename.empty() || filename == "auto")
        {
            throw std::runtime_error("Sternheimer SIAB production provenance requires an explicit " + label
                                     + " filename, not an empty/auto placeholder.");
        }
        const std::string path = join_path(directory, filename);
        std::ifstream input(path.c_str(), std::ios::binary);
        if (!input)
        {
            throw std::runtime_error("Cannot open required Sternheimer SIAB " + label + " file: " + path);
        }
        paths.push_back(path);
    }
    return paths;
}

double require_single_primitive_rcut(const std::vector<double>& rcuts)
{
    if (rcuts.size() != 1 || !std::isfinite(rcuts.front()) || rcuts.front() <= 0.0)
    {
        throw std::invalid_argument(
            "out_sternheimer_siab requires exactly one explicit positive finite bessel_nao_rcut.");
    }
    return rcuts.front();
}

std::string resolve_executable_path()
{
#if defined(__linux__)
    std::vector<char> path(4096, '\0');
    const ssize_t size = ::readlink("/proc/self/exe", path.data(), path.size() - 1);
    if (size <= 0 || static_cast<std::size_t>(size) >= path.size())
    {
        throw std::runtime_error("Cannot resolve /proc/self/exe for required Sternheimer SIAB executable provenance.");
    }
    path[static_cast<std::size_t>(size)] = '\0';
    return std::string(path.data());
#else
    throw std::runtime_error("Required Sternheimer SIAB executable provenance currently needs Linux /proc/self/exe.");
#endif
}

std::string require_source_commit(const std::string& compiled_commit)
{
    const std::string selected = normalized_commit_token(compiled_commit);
    if (!is_hex_string(selected, 40) && !is_hex_string(selected, 64))
    {
        throw std::runtime_error(
            "out_sternheimer_siab requires a build configured from an exact Git checkout with full COMMIT metadata.");
    }
    std::string normalized = selected;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return normalized;
}

} // namespace sternheimer_siab
} // namespace module_ri
