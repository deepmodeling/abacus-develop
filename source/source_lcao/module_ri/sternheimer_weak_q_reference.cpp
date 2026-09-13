#include "source_lcao/module_ri/sternheimer_weak_q_reference.h"
#include "source_lcao/module_ri/sternheimer_siab_provenance.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace ModuleRI
{
namespace
{
using Sha256 = module_ri::sternheimer_siab::Sha256;
constexpr unsigned char magic[8] = {'S', 'T', 'W', 'Q', 'R', 'E', 'F', 0};
constexpr std::uint32_t schema = 1;
constexpr std::uint64_t checksum_size = 64;
constexpr std::uint64_t record_header_size = 94;
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
              "Weak-q references require IEEE binary64");
static_assert(sizeof(int) == 4, "Weak-q record indices require 32-bit int");

[[noreturn]] void fail(const std::string& message)
{
    throw std::runtime_error("Sternheimer weak-q reference: " + message);
}

[[noreturn]] void io_fail(const std::string& operation)
{
    const int error = errno;
    fail(operation + ": " + std::strerror(error));
}

struct Budget
{
    std::uint64_t remaining;

    void charge(const std::uint64_t count, const std::uint64_t width)
    {
        if (width != 0 && count > remaining / width) fail("size/allocation limit exceeded");
        remaining -= count * width;
    }
};

void validate_arguments(const std::string& path, const std::string& contract,
                        const SternheimerWeakQReferenceDimensions& d,
                        const SternheimerWeakQReferenceLimits& limits)
{
    if (path.empty() || path.find('\0') != std::string::npos) fail("invalid path");
    if (contract.empty() || contract.size() > 65536) fail("invalid physical contract hash");
    if (d.fine_grid_size == 0 || d.full_record_count == 0 || d.orbital_count == 0)
        fail("dimensions must be positive");
    if (limits.max_file_bytes < checksum_size || limits.max_allocation_bytes == 0
        || limits.max_bands_per_record == 0)
        fail("invalid limits");
    for (const auto count : {d.fine_grid_size, d.full_record_count, d.orbital_count})
        if (count > std::numeric_limits<std::size_t>::max()) fail("dimension exceeds address space");
}

void charge_record(Budget& file, Budget& allocation, const std::uint64_t occupied,
                   const std::uint64_t unoccupied, const std::uint64_t orbitals,
                   const SternheimerWeakQReferenceLimits& limits)
{
    if (occupied > limits.max_bands_per_record || unoccupied > limits.max_bands_per_record - occupied)
        fail("band-count limit exceeded");
    file.charge(1, record_header_size);
    file.charge(occupied, 16); // Eigenvalues and occupations.
    file.charge(unoccupied, 8);
    allocation.charge(occupied, 2 * sizeof(double));
    allocation.charge(unoccupied, sizeof(double));
    const auto bands = occupied + unoccupied; // Addition bounded by the band limit above.
    if (bands != 0 && orbitals > std::numeric_limits<std::uint64_t>::max() / bands)
        fail("coefficient dimension overflow");
    const auto elements = bands * orbitals;
    file.charge(elements, 16);
    allocation.charge(bands, sizeof(std::vector<std::complex<double>>));
    allocation.charge(elements, sizeof(std::complex<double>));
}

class File
{
  public:
    std::FILE* stream = nullptr;
    std::string temporary;
    File() = default;
    File(const File&) = delete;
    File& operator=(const File&) = delete;
    ~File()
    {
        if (stream) std::fclose(stream);
        if (!temporary.empty()) ::unlink(temporary.c_str());
    }

    void attach(const int fd, const char* mode)
    {
        stream = ::fdopen(fd, mode);
        if (!stream)
        {
            const int error = errno;
            ::close(fd);
            errno = error;
            io_fail("fdopen");
        }
    }
};

class Writer
{
  public:
    Writer(std::FILE* stream, const std::uint64_t limit) : stream_(stream), budget_{limit} {}

    void bytes(const unsigned char* data, const std::size_t size)
    {
        budget_.charge(size, 1);
        if (std::fwrite(data, 1, size, stream_) != size) io_fail("write");
        sha_.update(data, size);
    }

    void integer(std::uint64_t value, const unsigned width = 8)
    {
        unsigned char data[8];
        for (unsigned i = 0; i < width; ++i)
        {
            data[i] = value & 255;
            value >>= 8;
        }
        bytes(data, width);
    }

    void real(const double value)
    {
        if (!std::isfinite(value)) fail("nonfinite value");
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        integer(bits);
    }

    void reals(const std::vector<double>& values)
    {
        for (const double value : values) real(value);
    }

    void coefficients(const std::vector<std::vector<std::complex<double>>>& values)
    {
        for (const auto& row : values)
            for (const auto& value : row)
            {
                real(value.real());
                real(value.imag());
            }
    }

    void seal()
    {
        const auto digest = sha_.finish();
        budget_.charge(digest.size(), 1);
        if (std::fwrite(digest.data(), 1, digest.size(), stream_) != digest.size()) io_fail("write checksum");
    }

  private:
    std::FILE* stream_;
    Budget budget_;
    Sha256 sha_;
};

class Reader
{
  public:
    Reader(std::FILE* stream, const std::uint64_t payload_size) : stream_(stream), remaining_(payload_size) {}

    std::uint64_t remaining() const { return remaining_; }

    void bytes(unsigned char* data, const std::size_t size)
    {
        if (size > remaining_) fail("truncated payload");
        if (std::fread(data, 1, size, stream_) != size) fail("truncated payload or read error");
        remaining_ -= size;
        sha_.update(data, size);
    }

    std::uint64_t integer(const unsigned width = 8)
    {
        unsigned char data[8];
        bytes(data, width);
        std::uint64_t value = 0;
        for (unsigned i = 0; i < width; ++i) value |= std::uint64_t(data[i]) << (8 * i);
        return value;
    }

    int index()
    {
        const auto value = integer(4);
        // Avoid implementation-defined unsigned-to-signed conversion.
        return value <= INT32_MAX ? static_cast<int>(value)
                                 : static_cast<int>(-1 - static_cast<std::int64_t>(UINT32_MAX - value));
    }

    bool boolean()
    {
        const auto value = integer(1);
        if (value > 1) fail("invalid boolean");
        return value == 1;
    }

    double real()
    {
        const auto bits = integer();
        double value;
        std::memcpy(&value, &bits, sizeof(value));
        if (!std::isfinite(value)) fail("nonfinite value");
        return value;
    }

    void reals(std::vector<double>& values, const std::uint64_t count)
    {
        if (count > remaining_ / 8 || count > values.max_size()) fail("invalid real-vector size");
        values.resize(static_cast<std::size_t>(count));
        for (auto& value : values) value = real();
    }

    void coefficients(std::vector<std::vector<std::complex<double>>>& values,
                      const std::uint64_t rows, const std::uint64_t columns)
    {
        if (rows > values.max_size()) fail("invalid coefficient row count");
        values.resize(static_cast<std::size_t>(rows));
        for (auto& row : values)
        {
            if (columns > remaining_ / 16 || columns > row.max_size()) fail("invalid coefficient width");
            row.resize(static_cast<std::size_t>(columns));
            for (auto& value : row)
            {
                const double re = real();
                const double im = real();
                value = {re, im};
            }
        }
    }

    void verify()
    {
        if (remaining_ != 0) fail("trailing payload bytes");
        const auto expected = sha_.finish();
        char actual[checksum_size];
        if (std::fread(actual, 1, sizeof(actual), stream_) != sizeof(actual)) fail("truncated checksum");
        if (std::memcmp(actual, expected.data(), sizeof(actual)) != 0) fail("SHA256 mismatch");
        if (std::fgetc(stream_) != EOF || std::ferror(stream_)) fail("trailing bytes or read error");
    }

  private:
    std::FILE* stream_;
    std::uint64_t remaining_;
    Sha256 sha_;
};

void validate_record_shape(const SternheimerLCAOOccupiedKPoint& r, const std::uint64_t orbitals)
{
    if (r.eigenvalues.size() != r.occupations.size() || r.eigenvalues.size() != r.coefficients.size()
        || r.unoccupied_eigenvalues.size() != r.unoccupied_coefficients.size())
        fail("inconsistent band dimensions");
    for (const auto* matrix : {&r.coefficients, &r.unoccupied_coefficients})
        for (const auto& row : *matrix)
            if (row.size() != orbitals) fail("inconsistent orbital dimensions");
}

void write_record(Writer& out, const SternheimerLCAOOccupiedKPoint& r)
{
    for (const int value : {r.local_k_index, r.global_k_index, r.zero_order_k_index,
                           r.symmetry_spatial_isym, r.spin_index})
        out.integer(static_cast<std::uint32_t>(value), 4);
    out.integer(r.symmetry_time_reversal ? 1 : 0, 1);
    out.integer(r.has_grid_kpoint_override ? 1 : 0, 1);
    for (const double value : r.kpoint) out.real(value);
    for (const double value : r.grid_kpoint) out.real(value);
    out.real(r.kweight);
    out.integer(r.eigenvalues.size());
    out.integer(r.unoccupied_eigenvalues.size());
    out.reals(r.eigenvalues);
    out.reals(r.occupations);
    out.coefficients(r.coefficients);
    out.reals(r.unoccupied_eigenvalues);
    out.coefficients(r.unoccupied_coefficients);
}

void read_record(Reader& in, SternheimerLCAOOccupiedKPoint& r, Budget& allocation,
                 const std::uint64_t orbitals, const SternheimerWeakQReferenceLimits& limits)
{
    const auto available = in.remaining();
    r.local_k_index = in.index();
    r.global_k_index = in.index();
    r.zero_order_k_index = in.index();
    r.symmetry_spatial_isym = in.index();
    r.spin_index = in.index();
    r.symmetry_time_reversal = in.boolean();
    r.has_grid_kpoint_override = in.boolean();
    for (auto& value : r.kpoint) value = in.real();
    for (auto& value : r.grid_kpoint) value = in.real();
    r.kweight = in.real();
    const auto occupied = in.integer();
    const auto unoccupied = in.integer();
    // Check encoded size and cumulative decoded storage BEFORE any band allocation.
    Budget file{available};
    charge_record(file, allocation, occupied, unoccupied, orbitals, limits);
    in.reals(r.eigenvalues, occupied);
    in.reals(r.occupations, occupied);
    in.coefficients(r.coefficients, occupied, orbitals);
    in.reals(r.unoccupied_eigenvalues, unoccupied);
    in.coefficients(r.unoccupied_coefficients, unoccupied, orbitals);
}
} // namespace

void write_sternheimer_weak_q_reference(
    const std::string& path, const std::string& contract_hash,
    const SternheimerWeakQReferenceDimensions& d, const std::vector<double>& fine_potential,
    const std::vector<SternheimerLCAOOccupiedKPoint>& full_records,
    const SternheimerWeakQReferenceLimits& limits)
{
    validate_arguments(path, contract_hash, d, limits);
    if (fine_potential.size() != d.fine_grid_size || full_records.size() != d.full_record_count)
        fail("potential/record dimensions mismatch");
    Budget file_budget{limits.max_file_bytes};
    Budget allocation{limits.max_allocation_bytes};
    file_budget.charge(1, 8 + 4 + 8 + 24 + checksum_size);
    file_budget.charge(contract_hash.size(), 1);
    file_budget.charge(d.fine_grid_size, 8);
    allocation.charge(d.fine_grid_size, sizeof(double));
    allocation.charge(d.full_record_count, sizeof(SternheimerLCAOOccupiedKPoint));
    for (const auto& r : full_records)
    {
        validate_record_shape(r, d.orbital_count);
        charge_record(file_budget, allocation, r.eigenvalues.size(), r.unoccupied_eigenvalues.size(),
                      d.orbital_count, limits);
    }

    File file;
    // Allocate the cleanup path before creating the file; mkstemp mutates only its suffix.
    file.temporary = path + ".tmp.XXXXXX";
    const int fd = ::mkstemp(&file.temporary[0]);
    if (fd < 0)
    {
        file.temporary.clear();
        io_fail("create temporary file");
    }
    file.attach(fd, "wb");
    Writer out(file.stream, limits.max_file_bytes);
    out.bytes(magic, sizeof(magic));
    out.integer(schema, 4);
    out.integer(contract_hash.size());
    out.bytes(reinterpret_cast<const unsigned char*>(contract_hash.data()), contract_hash.size());
    out.integer(d.fine_grid_size);
    out.integer(d.full_record_count);
    out.integer(d.orbital_count);
    out.reals(fine_potential);
    for (const auto& r : full_records) write_record(out, r);
    out.seal();
    if (std::fflush(file.stream) != 0) io_fail("flush temporary file");
    if (::fsync(::fileno(file.stream)) != 0) io_fail("sync temporary file");
    const int close_result = std::fclose(file.stream);
    file.stream = nullptr;
    if (close_result != 0) io_fail("close temporary file");
    // link(), unlike rename(), has atomic no-replace semantics on POSIX.
    if (::link(file.temporary.c_str(), path.c_str()) != 0) io_fail("publish new file");
}

SternheimerWeakQReference read_sternheimer_weak_q_reference(
    const std::string& path, const std::string& contract_hash,
    const SternheimerWeakQReferenceDimensions& d, const SternheimerWeakQReferenceLimits& limits)
{
    validate_arguments(path, contract_hash, d, limits);
    File file;
    // Nonblocking open avoids hanging on a FIFO before the regular-file check.
    const int fd = ::open(path.c_str(), O_RDONLY | O_NONBLOCK);
    if (fd < 0) io_fail("open reference");
    file.attach(fd, "rb");
    struct stat st;
    if (::fstat(fd, &st) != 0) io_fail("stat reference");
    if (!S_ISREG(st.st_mode) || st.st_size < static_cast<off_t>(checksum_size)) fail("invalid reference file");
    const auto size = static_cast<std::uint64_t>(st.st_size);
    if (size > limits.max_file_bytes) fail("file-size limit exceeded");
    Reader in(file.stream, size - checksum_size);
    unsigned char found_magic[sizeof(magic)];
    in.bytes(found_magic, sizeof(found_magic));
    if (std::memcmp(found_magic, magic, sizeof(magic)) != 0) fail("invalid magic");
    if (in.integer(4) != schema) fail("unsupported schema");
    if (in.integer() != contract_hash.size()) fail("physical contract mismatch");
    unsigned char chunk[4096];
    for (std::size_t offset = 0; offset < contract_hash.size();)
    {
        const auto length = std::min(sizeof(chunk), contract_hash.size() - offset);
        in.bytes(chunk, length);
        if (std::memcmp(chunk, contract_hash.data() + offset, length) != 0) fail("physical contract mismatch");
        offset += length;
    }
    const auto grid = in.integer();
    const auto count = in.integer();
    const auto orbitals = in.integer();
    if (grid != d.fine_grid_size || count != d.full_record_count || orbitals != d.orbital_count)
        fail("physical dimensions mismatch");

    Budget encoded{in.remaining()};
    encoded.charge(grid, 8);
    encoded.charge(count, record_header_size);
    Budget allocation{limits.max_allocation_bytes};
    allocation.charge(grid, sizeof(double));
    allocation.charge(count, sizeof(SternheimerLCAOOccupiedKPoint));
    SternheimerWeakQReference result;
    if (count > result.full_records.max_size()) fail("invalid record count");
    in.reals(result.fine_potential, grid);
    result.full_records.resize(static_cast<std::size_t>(count));
    for (auto& r : result.full_records) read_record(in, r, allocation, orbitals, limits);
    in.verify();
    return result;
}

} // namespace ModuleRI
