//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdint>
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
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include "rpa_lri.h"
#include "exx_lri.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_base/global_function.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/utils/spectrum_mo.hpp"

#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace RpaLriDetail
{
constexpr int LIBRPA_COULOMB_V1_MARKER = -20129433;
constexpr int LIBRPA_LRICOEF_V1_MARKER = -10267453;
constexpr int LIBRPA_SHRINK_SINVS_V1_MARKER = -30241621;
constexpr int LIBRPA_KS_EIGENVECTOR_V1_MARKER = -12345679;
constexpr int LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE = 28;
constexpr int LIBRPA_COULOMB_V1_COMPLEX_FLAG = 1;

static_assert(sizeof(std::complex<double>) == 2 * sizeof(double),
              "LibRPA v1 binary output expects complex<double> as two doubles.");

inline void trim_malloc_cache()
{
#if defined(__GLIBC__)
    malloc_trim(0);
#endif
}
inline bool debug_dump_exx_ao_enabled()
{
    const char* env = std::getenv("ABACUS_DUMP_EXX_AO");
    if (env == nullptr)
    {
        return false;
    }
    const std::string value(env);
    return !(value.empty() || value == "0" || value == "f" || value == "F"
             || value == "false" || value == "FALSE");
}

inline std::size_t coulomb_atom_pair_index(const std::size_t I, const std::size_t J, const std::size_t natoms)
{
    if (I > J)
    {
        throw std::runtime_error("LibRPA v1 Coulomb output expects upper-triangular atom pairs.");
    }
    return I * natoms - I * (I - 1) / 2 + (J - I);
}

inline void checked_write(std::ofstream& ofs, const void* data, const std::size_t bytes, const std::string& filename)
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
inline void write_scalar(std::ofstream& ofs, const Value& value, const std::string& filename)
{
    checked_write(ofs, &value, sizeof(Value), filename);
}

inline unsigned long long checked_mul_u64(const unsigned long long lhs,
                                          const unsigned long long rhs,
                                          const std::string& context)
{
    if (lhs != 0 && rhs > std::numeric_limits<unsigned long long>::max() / lhs)
    {
        throw std::runtime_error(context + " exceeds uint64_t range.");
    }
    return lhs * rhs;
}

inline unsigned long long checked_add_u64(const unsigned long long lhs,
                                          const unsigned long long rhs,
                                          const std::string& context)
{
    if (rhs > std::numeric_limits<unsigned long long>::max() - lhs)
    {
        throw std::runtime_error(context + " exceeds uint64_t range.");
    }
    return lhs + rhs;
}

inline std::int64_t checked_i64_from_u64(const unsigned long long value, const std::string& context)
{
    if (value > static_cast<unsigned long long>(std::numeric_limits<std::int64_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int64_t range.");
    }
    return static_cast<std::int64_t>(value);
}

inline std::int32_t checked_i32_from_size(const std::size_t value, const std::string& context)
{
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    {
        throw std::runtime_error(context + " exceeds int32_t range.");
    }
    return static_cast<std::int32_t>(value);
}

inline std::int32_t checked_i32_from_int(const int value, const std::string& context)
{
    if (value < 0)
    {
        throw std::runtime_error(context + " is negative.");
    }
    return checked_i32_from_size(static_cast<std::size_t>(value), context);
}

#ifdef __MPI
inline std::string mpi_error_string(const int error_code)
{
    char error_buffer[MPI_MAX_ERROR_STRING] = {};
    int error_length = 0;
    if (MPI_Error_string(error_code, error_buffer, &error_length) != MPI_SUCCESS)
    {
        return "MPI error " + std::to_string(error_code);
    }
    return std::string(error_buffer, static_cast<std::size_t>(error_length));
}

inline void collective_mpi_check(const MPI_Comm mpi_comm,
                                 const int local_error,
                                 const std::string& context)
{
    int mpi_rank = 0;
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    const int no_failure = std::numeric_limits<int>::max();
    const int local_failed_rank = local_error == MPI_SUCCESS ? no_failure : mpi_rank;
    int first_failed_rank = no_failure;
    const int reduce_error
        = MPI_Allreduce(&local_failed_rank, &first_failed_rank, 1, MPI_INT, MPI_MIN, mpi_comm);
    if (reduce_error != MPI_SUCCESS)
    {
        throw std::runtime_error(context + ": failed to synchronize MPI error status: "
                                 + mpi_error_string(reduce_error));
    }
    if (first_failed_rank == no_failure)
    {
        return;
    }

    int shared_error = mpi_rank == first_failed_rank ? local_error : MPI_SUCCESS;
    const int broadcast_error = MPI_Bcast(&shared_error, 1, MPI_INT, first_failed_rank, mpi_comm);
    if (broadcast_error != MPI_SUCCESS)
    {
        throw std::runtime_error(context + ": failed to broadcast MPI error status: "
                                 + mpi_error_string(broadcast_error));
    }
    throw std::runtime_error(context + ": " + mpi_error_string(shared_error));
}

inline void collective_require(const MPI_Comm mpi_comm,
                               const bool local_condition,
                               const std::string& context)
{
    int mpi_rank = 0;
    MPI_Comm_rank(mpi_comm, &mpi_rank);
    const int no_failure = std::numeric_limits<int>::max();
    const int local_failed_rank = local_condition ? no_failure : mpi_rank;
    int first_failed_rank = no_failure;
    const int reduce_error
        = MPI_Allreduce(&local_failed_rank, &first_failed_rank, 1, MPI_INT, MPI_MIN, mpi_comm);
    if (reduce_error != MPI_SUCCESS)
    {
        throw std::runtime_error(context + ": failed to synchronize validation status: "
                                 + mpi_error_string(reduce_error));
    }
    if (first_failed_rank != no_failure)
    {
        throw std::runtime_error(context + " (first failing MPI rank "
                                 + std::to_string(first_failed_rank) + ").");
    }
}

inline MPI_Aint checked_mpi_aint_from_u64(const unsigned long long value, const std::string& context)
{
    if (value > static_cast<unsigned long long>(std::numeric_limits<MPI_Aint>::max()))
    {
        throw std::runtime_error(context + " exceeds MPI_Aint range.");
    }
    return static_cast<MPI_Aint>(value);
}

inline MPI_Offset checked_mpi_offset_from_u64(const unsigned long long value, const std::string& context)
{
    if (value > static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max()))
    {
        throw std::runtime_error(context + " exceeds MPI_Offset range.");
    }
    return static_cast<MPI_Offset>(value);
}

struct KSEigenvectorMpiLayout
{
    MPI_Datatype filetype = MPI_C_DOUBLE_COMPLEX;
    bool free_filetype = false;
    unsigned long long local_count = 0;
    unsigned long long max_local_count = 0;
};

inline KSEigenvectorMpiLayout make_ks_eigenvector_mpi_layout(const MPI_Comm mpi_comm,
                                                              const Parallel_Orbitals& parav,
                                                              const int nbands,
                                                              const int nbasis_wfc,
                                                              const bool is_soc,
                                                              const int spinor_component)
{
    KSEigenvectorMpiLayout layout;
    std::vector<int> block_lengths;
    std::vector<MPI_Aint> block_displacements;
    bool local_valid = nbands >= 0 && nbasis_wfc >= 0 && parav.ncol_bands >= 0
                       && parav.ncol_bands <= parav.get_col_size()
                       && spinor_component >= 0 && spinor_component < (is_soc ? 2 : 1);
    const int spatial_basis = is_soc ? nbasis_wfc / 2 : nbasis_wfc;

    try
    {
        bool have_run = false;
        unsigned long long run_first_index = 0;
        unsigned long long previous_index = 0;
        int run_length = 0;

        const auto flush_run = [&]()
        {
            if (!have_run)
            {
                return;
            }
            const unsigned long long byte_displacement
                = checked_mul_u64(run_first_index,
                                  static_cast<unsigned long long>(sizeof(std::complex<double>)),
                                  "KS eigenvector MPI file displacement");
            block_lengths.push_back(run_length);
            block_displacements.push_back(
                checked_mpi_aint_from_u64(byte_displacement, "KS eigenvector MPI file displacement"));
            have_run = false;
            run_length = 0;
        };

        for (int ib_local = 0; local_valid && ib_local < parav.ncol_bands; ++ib_local)
        {
            const int ib_global = parav.local2global_col(ib_local);
            for (int ir_local = 0; ir_local < parav.get_row_size(); ++ir_local)
            {
                const int iw_global = parav.local2global_row(ir_local);
                if (is_soc && iw_global % 2 != spinor_component)
                {
                    continue;
                }

                const int iw_file = is_soc ? iw_global / 2 : iw_global;
                const unsigned long long band_offset
                    = checked_mul_u64(static_cast<unsigned long long>(ib_global),
                                      static_cast<unsigned long long>(spatial_basis),
                                      "KS eigenvector MPI file index");
                const unsigned long long file_index
                    = checked_add_u64(band_offset,
                                      static_cast<unsigned long long>(iw_file),
                                      "KS eigenvector MPI file index");
                if (!have_run)
                {
                    have_run = true;
                    run_first_index = file_index;
                    previous_index = file_index;
                    run_length = 1;
                }
                else if (file_index == previous_index + 1
                         && run_length < std::numeric_limits<int>::max())
                {
                    ++run_length;
                    previous_index = file_index;
                }
                else
                {
                    flush_run();
                    have_run = true;
                    run_first_index = file_index;
                    previous_index = file_index;
                    run_length = 1;
                }
                layout.local_count = checked_add_u64(layout.local_count,
                                                     1,
                                                     "Local KS eigenvector MPI element count");
            }
        }
        flush_run();
    }
    catch (const std::exception&)
    {
        local_valid = false;
    }

    local_valid = local_valid
                  && block_lengths.size()
                         <= static_cast<std::size_t>(std::numeric_limits<int>::max());
    collective_require(mpi_comm, local_valid, "Invalid local KS eigenvector MPI layout");

    unsigned long long global_count = 0;
    collective_mpi_check(
        mpi_comm,
        MPI_Allreduce(&layout.local_count,
                      &global_count,
                      1,
                      MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM,
                      mpi_comm),
        "Failed to validate KS eigenvector MPI ownership");
    const unsigned long long expected_count
        = checked_mul_u64(static_cast<unsigned long long>(nbands),
                          static_cast<unsigned long long>(spatial_basis),
                          "Global KS eigenvector MPI element count");
    collective_require(mpi_comm,
                       global_count == expected_count,
                       "KS eigenvector MPI ownership does not cover the global matrix exactly");

    collective_mpi_check(
        mpi_comm,
        MPI_Allreduce(&layout.local_count,
                      &layout.max_local_count,
                      1,
                      MPI_UNSIGNED_LONG_LONG,
                      MPI_MAX,
                      mpi_comm),
        "Failed to determine the KS eigenvector MPI chunk count");

    int type_create_error = MPI_SUCCESS;
    if (layout.local_count > 0)
    {
        type_create_error = MPI_Type_create_hindexed(static_cast<int>(block_lengths.size()),
                                                     block_lengths.data(),
                                                     block_displacements.data(),
                                                     MPI_C_DOUBLE_COMPLEX,
                                                     &layout.filetype);
    }
    collective_mpi_check(mpi_comm,
                         type_create_error,
                         "Failed to create KS eigenvector MPI file type");

    int type_commit_error = MPI_SUCCESS;
    if (layout.local_count > 0)
    {
        type_commit_error = MPI_Type_commit(&layout.filetype);
    }
    collective_mpi_check(mpi_comm,
                         type_commit_error,
                         "Failed to commit KS eigenvector MPI file type");
    layout.free_filetype = layout.local_count > 0;
    return layout;
}

inline int free_ks_eigenvector_mpi_layout(KSEigenvectorMpiLayout& layout)
{
    if (!layout.free_filetype)
    {
        return MPI_SUCCESS;
    }
    layout.free_filetype = false;
    return MPI_Type_free(&layout.filetype);
}

template <typename Value>
inline void append_binary_scalar(std::vector<char>& bytes, const Value& value)
{
    const std::size_t old_size = bytes.size();
    bytes.resize(old_size + sizeof(Value));
    std::memcpy(bytes.data() + old_size, &value, sizeof(Value));
}

struct KSEigenvectorV1Metadata
{
    std::vector<char> bytes;
    std::vector<unsigned long long> payload_offsets;
    unsigned long long component_bytes = 0;
    unsigned long long total_bytes = 0;
};

inline KSEigenvectorV1Metadata make_ks_eigenvector_v1_metadata(const int nks_tot,
                                                               const int nspins,
                                                               const int ncomponents,
                                                               const int nbands,
                                                               const int nbasis_wfc,
                                                               const int spatial_basis)
{
    KSEigenvectorV1Metadata metadata;
    const unsigned long long header_bytes
        = checked_mul_u64(6,
                          static_cast<unsigned long long>(sizeof(std::int32_t)),
                          "KS eigenvector v1 header size");
    const unsigned long long record_bytes
        = checked_add_u64(static_cast<unsigned long long>(sizeof(std::int32_t)),
                          static_cast<unsigned long long>(sizeof(std::int64_t)),
                          "KS eigenvector v1 directory record size");
    const unsigned long long directory_bytes
        = checked_mul_u64(static_cast<unsigned long long>(nks_tot),
                          record_bytes,
                          "KS eigenvector v1 directory size");
    const unsigned long long payload_begin
        = checked_add_u64(header_bytes, directory_bytes, "KS eigenvector v1 metadata size");
    const unsigned long long component_values
        = checked_mul_u64(static_cast<unsigned long long>(nbands),
                          static_cast<unsigned long long>(spatial_basis),
                          "KS eigenvector v1 component size");
    metadata.component_bytes
        = checked_mul_u64(component_values,
                          static_cast<unsigned long long>(sizeof(std::complex<double>)),
                          "KS eigenvector v1 component byte size");
    const unsigned long long kpoint_bytes
        = checked_mul_u64(static_cast<unsigned long long>(ncomponents),
                          metadata.component_bytes,
                          "KS eigenvector v1 k-point byte size");
#endif
