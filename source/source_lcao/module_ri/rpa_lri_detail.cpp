#include "rpa_lri_detail.h"

#include <RI/global/Tensor.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace RpaLriDetail
{
void trim_malloc_cache()
{
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}
bool debug_dump_exx_ao_enabled()
{
    const char* env = std::getenv("ABACUS_DUMP_EXX_AO");
    if (env == nullptr)
    {
        return false;
    }
    const std::string value(env);
    return !(value.empty() || value == "0" || value == "f" || value == "F" || value == "false" || value == "FALSE");
}

std::size_t coulomb_atom_pair_index(const std::size_t I, const std::size_t J, const std::size_t natoms)
{
    if (I > J)
    {
        throw std::runtime_error("LibRPA v1 Coulomb output expects upper-triangular atom pairs.");
    }
    return I * natoms - I * (I - 1) / 2 + (J - I);
}

void checked_write(std::ofstream& ofs, const void* data, const std::size_t bytes, const std::string& filename)
{
    const char* ptr = reinterpret_cast<const char*>(data);
    std::size_t bytes_left = bytes;
    const std::size_t max_chunk = static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max());
    while (bytes_left > 0)
    {
        const std::size_t chunk = std::min(bytes_left, max_chunk);
        ofs.write(ptr, static_cast<std::streamsize>(chunk));
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to write " + filename);
        }
        ptr += chunk;
        bytes_left -= chunk;
    }
}

template <typename Value>
void write_scalar(std::ofstream& ofs, const Value& value, const std::string& filename)
{
    checked_write(ofs, &value, sizeof(Value), filename);
}

unsigned long long checked_mul_u64(const unsigned long long lhs,
                                   const unsigned long long rhs,
                                   const std::string& context)
{
    if (lhs != 0 && rhs > std::numeric_limits<unsigned long long>::max() / lhs)
    {
        throw std::runtime_error(context + " exceeds uint64_t range.");
    }
    return lhs * rhs;
}

unsigned long long checked_add_u64(const unsigned long long lhs,
                                   const unsigned long long rhs,
                                   const std::string& context)
{
    if (rhs > std::numeric_limits<unsigned long long>::max() - lhs)
    {
        throw std::runtime_error(context + " exceeds uint64_t range.");
    }
    return lhs + rhs;
}

std::int64_t checked_i64_from_u64(const unsigned long long value, const std::string& context)
{
    if (value > static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int64_t range.");
    }
    return static_cast<std::int64_t>(value);
}

std::int32_t checked_i32_from_size(const std::size_t value, const std::string& context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int32_t range.");
    }
    return static_cast<std::int32_t>(value);
}

std::int32_t checked_i32_from_int(const int value, const std::string& context)
{
    if (value < 0)
    {
        throw std::runtime_error(context + " is negative.");
    }
    return checked_i32_from_size(static_cast<std::size_t>(value), context);
}
int sum_int_vector(const std::vector<int>& values)
{
    int sum = 0;
    for (const int value: values)
    {
        if (value > std::numeric_limits<int>::max() - sum)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        sum += value;
    }
    return sum;
}
double real_as_double(const std::complex<double>& value)
{
    return value.real();
}
double real_as_double(double value)
{
    return value;
}

bool has_valid_matrix_shape(const RI::Tensor<double>& tensor)
{
    return tensor.shape.size() == 2 && tensor.shape[0] > 0 && tensor.shape[1] > 0;
}

void abort_output(MPI_Comm comm, const std::string& message)
{
    std::cerr << "RPA producer output failed: " << message << std::endl;
#ifdef __MPI
    MPI_Abort(comm, 1);
#endif
    std::exit(1);
}
template void write_scalar<std::int32_t>(std::ofstream&, const std::int32_t&, const std::string&);
template void write_scalar<std::int64_t>(std::ofstream&, const std::int64_t&, const std::string&);
template void write_scalar<double>(std::ofstream&, const double&, const std::string&);
} // namespace RpaLriDetail
