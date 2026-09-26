//=======================
// AUTHOR : Rong Shi
// DATE :   2022-12-09
//=======================

#ifndef RPA_LRI_HPP
#define RPA_LRI_HPP
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"
#include "source_cell/module_symmetry/symm_rot_spin.h"

#include "rpa_lri.h"
#include "librpa_2d_coulomb_head.h"
#include "librpa_bz_sampling.h"
#include "librpa_stru_symmetry.h"
#include "librpa_stru_units.h"
#include "rpa_abfs_preorthogonalization.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_lcao/module_lr/utils/spectrum_mo.hpp"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_restart/restart_exx_csr.h"

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
constexpr int LIBRPA_ABF_OVERLAP_V1_MARKER = -40817329;
constexpr int LIBRPA_ABF_OVERLAP_V1_VERSION = 1;
constexpr int LIBRPA_ABF_OVERLAP_V1_KIND_ACTIVE = 1;

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

inline bool ewald_component_output_enabled()
{
    const char* env = std::getenv("ABACUS_RPA_EWALD_COMPONENTS");
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
    const unsigned long long payload_bytes
        = checked_mul_u64(static_cast<unsigned long long>(nks_tot),
                          kpoint_bytes,
                          "KS eigenvector v1 payload byte size");
    metadata.total_bytes
        = checked_add_u64(payload_begin, payload_bytes, "KS eigenvector v1 file size");

    const std::int32_t marker = LIBRPA_KS_EIGENVECTOR_V1_MARKER;
    const std::int32_t kind = LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE;
    const std::int32_t nkpoints_i32 = checked_i32_from_int(nks_tot, "KS eigenvector k-point count");
    const std::int32_t nspins_i32 = checked_i32_from_int(nspins, "KS eigenvector spin count");
    const std::int32_t nstates_i32 = checked_i32_from_int(nbands, "KS eigenvector state count");
    const std::int32_t nbasis_i32 = checked_i32_from_int(nbasis_wfc, "KS eigenvector basis count");
    append_binary_scalar(metadata.bytes, marker);
    append_binary_scalar(metadata.bytes, kind);
    append_binary_scalar(metadata.bytes, nkpoints_i32);
    append_binary_scalar(metadata.bytes, nspins_i32);
    append_binary_scalar(metadata.bytes, nstates_i32);
    append_binary_scalar(metadata.bytes, nbasis_i32);

    metadata.payload_offsets.reserve(static_cast<std::size_t>(nks_tot));
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const unsigned long long kpoint_offset
            = checked_mul_u64(static_cast<unsigned long long>(ik),
                              kpoint_bytes,
                              "KS eigenvector v1 k-point offset");
        const unsigned long long payload_offset
            = checked_add_u64(payload_begin, kpoint_offset, "KS eigenvector v1 payload offset");
        metadata.payload_offsets.push_back(payload_offset);
        const std::int32_t ik_file = checked_i32_from_int(ik + 1, "KS eigenvector k-point index");
        const std::int64_t payload_offset_i64
            = checked_i64_from_u64(payload_offset, "KS eigenvector v1 payload offset");
        append_binary_scalar(metadata.bytes, ik_file);
        append_binary_scalar(metadata.bytes, payload_offset_i64);
    }
    if (metadata.bytes.size() != static_cast<std::size_t>(payload_begin))
    {
        throw std::runtime_error("KS eigenvector v1 metadata size is inconsistent.");
    }
    return metadata;
}

struct KSEigenvectorPackCursor
{
    int band_local = 0;
    int row_local = 0;
};

template <typename T>
inline int pack_ks_eigenvector_chunk(const Parallel_Orbitals& parav,
                                     const psi::Psi<T>& psi,
                                     const int psi_k,
                                     const bool is_soc,
                                     const int spinor_component,
                                     KSEigenvectorPackCursor& cursor,
                                     std::vector<std::complex<double>>& buffer,
                                     const int requested_count)
{
    int packed = 0;
    while (cursor.band_local < parav.ncol_bands && packed < requested_count)
    {
        while (cursor.row_local < psi.get_nbasis() && packed < requested_count)
        {
            const int ir_local = cursor.row_local++;
            const int iw_global = parav.local2global_row(ir_local);
            if (is_soc && iw_global % 2 != spinor_component)
            {
                continue;
            }
            buffer[static_cast<std::size_t>(packed)]
                = std::complex<double>(psi(psi_k, cursor.band_local, ir_local));
            ++packed;
        }
        if (cursor.row_local == psi.get_nbasis())
        {
            ++cursor.band_local;
            cursor.row_local = 0;
        }
    }
    return packed;
}

template <typename T>
inline void write_ks_eigenvector_v1_mpi(const MPI_Comm mpi_comm,
                                        const Parallel_Orbitals& parav,
                                        const psi::Psi<T>& psi,
                                        const int nks_tot,
                                        const int nspins,
                                        const int nspin_abacus,
                                        const int nbands,
                                        const int nbasis_wfc,
                                        const std::string& final_name)
{
    if (mpi_comm == MPI_COMM_NULL)
    {
        throw std::runtime_error("KS eigenvector MPI-IO writer received MPI_COMM_NULL.");
    }
    int mpi_rank = 0;
    collective_mpi_check(mpi_comm,
                         MPI_Comm_rank(mpi_comm, &mpi_rank),
                         "Failed to query the KS eigenvector MPI rank");

    const bool is_soc = nspin_abacus == 4;
    const int ncomponents = is_soc ? 2 : nspins;
    bool dimensions_valid
        = nks_tot >= 0 && nbands >= 0 && nbasis_wfc >= 0 && nspins > 0
          && (nspin_abacus == 1 || nspin_abacus == 2 || nspin_abacus == 4)
          && nspins == (nspin_abacus == 2 ? 2 : 1)
          && (!is_soc || nbasis_wfc % 2 == 0)
          && psi.get_nbasis() == parav.get_row_size()
          && psi.get_nbands() == parav.ncol_bands
          && parav.get_wfc_global_nbasis() == nbasis_wfc
          && parav.get_wfc_global_nbands() == nbands;
    if (dimensions_valid)
    {
        const unsigned long long required_psi_kpoints
            = checked_mul_u64(static_cast<unsigned long long>(nks_tot),
                              static_cast<unsigned long long>(nspins),
                              "KS eigenvector MPI source k-point count");
        dimensions_valid
            = required_psi_kpoints <= static_cast<unsigned long long>(psi.get_nk());
    }
    collective_require(mpi_comm, dimensions_valid, "Invalid KS eigenvector MPI-IO dimensions");

    const int spatial_basis = is_soc ? nbasis_wfc / 2 : nbasis_wfc;
    std::vector<KSEigenvectorMpiLayout> layouts;
    layouts.reserve(static_cast<std::size_t>(is_soc ? 2 : 1));
    for (int component = 0; component < (is_soc ? 2 : 1); ++component)
    {
        layouts.push_back(make_ks_eigenvector_mpi_layout(
            mpi_comm, parav, nbands, nbasis_wfc, is_soc, component));
    }

    const unsigned long long chunk_elements = 4ULL * 1024ULL * 1024ULL;
    unsigned long long local_buffer_elements = 0;
    for (const auto& layout: layouts)
    {
        local_buffer_elements = std::max(local_buffer_elements,
                                         std::min(layout.local_count, chunk_elements));
    }
    bool buffer_allocated = true;
    std::vector<std::complex<double>> buffer;
    try
    {
        buffer.resize(static_cast<std::size_t>(local_buffer_elements));
    }
    catch (const std::exception&)
    {
        buffer_allocated = false;
    }
    collective_require(mpi_comm,
                       buffer_allocated,
                       "Failed to allocate the bounded KS eigenvector MPI pack buffer");
    std::complex<double> dummy(0.0, 0.0);

    const KSEigenvectorV1Metadata metadata
        = make_ks_eigenvector_v1_metadata(
            nks_tot, nspins, ncomponents, nbands, nbasis_wfc, spatial_basis);
    collective_require(
        mpi_comm,
        metadata.bytes.size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "KS eigenvector MPI metadata exceeds the MPI count range");

    const std::string temporary_name = final_name + ".tmp";
    MPI_File file = MPI_FILE_NULL;
    bool file_opened = false;
    try
    {
        int mpi_error = MPI_File_open(mpi_comm,
                                      const_cast<char*>(temporary_name.c_str()),
                                      MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                      MPI_INFO_NULL,
                                      &file);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to open " + temporary_name);
        file_opened = true;
        collective_mpi_check(mpi_comm,
                             MPI_File_set_errhandler(file, MPI_ERRORS_RETURN),
                             "Failed to set the KS eigenvector MPI file error handler");
        collective_mpi_check(
            mpi_comm,
            MPI_File_set_size(
                file,
                checked_mpi_offset_from_u64(metadata.total_bytes, "KS eigenvector MPI file size")),
            "Failed to set the KS eigenvector MPI file size");

        const int metadata_count
            = mpi_rank == 0 ? static_cast<int>(metadata.bytes.size()) : 0;
        MPI_Status metadata_status;
        mpi_error = MPI_File_write_at_all(file,
                                          0,
                                          metadata_count == 0
                                              ? static_cast<void*>(&dummy)
                                              : static_cast<void*>(
                                                  const_cast<char*>(metadata.bytes.data())),
                                          metadata_count,
                                          MPI_BYTE,
                                          &metadata_status);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to write KS eigenvector MPI metadata");

        for (int ik = 0; ik < nks_tot; ++ik)
        {
            for (int component = 0; component < ncomponents; ++component)
            {
                KSEigenvectorMpiLayout& layout = layouts[is_soc ? component : 0];
                const unsigned long long component_offset
                    = checked_mul_u64(static_cast<unsigned long long>(component),
                                      metadata.component_bytes,
                                      "KS eigenvector MPI component offset");
                const unsigned long long displacement
                    = checked_add_u64(metadata.payload_offsets[static_cast<std::size_t>(ik)],
                                      component_offset,
                                      "KS eigenvector MPI view displacement");
                collective_mpi_check(
                    mpi_comm,
                    MPI_File_set_view(
                        file,
                        checked_mpi_offset_from_u64(displacement,
                                                    "KS eigenvector MPI view displacement"),
                        MPI_C_DOUBLE_COMPLEX,
                        layout.filetype,
                        const_cast<char*>("native"),
                        MPI_INFO_NULL),
                    "Failed to set the KS eigenvector MPI file view");

                KSEigenvectorPackCursor cursor;
                const int psi_k = is_soc ? ik : ik + nks_tot * component;
                unsigned long long chunk_begin = 0;
                while (chunk_begin < layout.max_local_count)
                {
                    const unsigned long long local_remaining
                        = chunk_begin < layout.local_count ? layout.local_count - chunk_begin : 0;
                    const int local_chunk_count = static_cast<int>(
                        std::min(local_remaining, chunk_elements));
                    const int packed = pack_ks_eigenvector_chunk(parav,
                                                                 psi,
                                                                 psi_k,
                                                                 is_soc,
                                                                 component,
                                                                 cursor,
                                                                 buffer,
                                                                 local_chunk_count);
                    collective_require(mpi_comm,
                                       packed == local_chunk_count,
                                       "Failed to pack the expected KS eigenvector MPI chunk");

                    MPI_Status payload_status;
                    mpi_error = MPI_File_write_at_all(
                        file,
                        checked_mpi_offset_from_u64(chunk_begin,
                                                    "KS eigenvector MPI chunk offset"),
                        local_chunk_count == 0 ? static_cast<void*>(&dummy)
                                               : static_cast<void*>(buffer.data()),
                        local_chunk_count,
                        MPI_C_DOUBLE_COMPLEX,
                        &payload_status);
                    collective_mpi_check(mpi_comm,
                                         mpi_error,
                                         "Failed to write a KS eigenvector MPI payload chunk");
                    if (layout.max_local_count - chunk_begin <= chunk_elements)
                    {
                        break;
                    }
                    chunk_begin = checked_add_u64(
                        chunk_begin, chunk_elements, "KS eigenvector MPI chunk offset");
                }
            }
        }

        mpi_error = MPI_File_sync(file);
        collective_mpi_check(mpi_comm, mpi_error, "Failed to sync the KS eigenvector MPI file");
        mpi_error = MPI_File_close(&file);
        file_opened = false;
        collective_mpi_check(mpi_comm, mpi_error, "Failed to close the KS eigenvector MPI file");

        int free_error = MPI_SUCCESS;
        for (auto& layout: layouts)
        {
            const int layout_error = free_ks_eigenvector_mpi_layout(layout);
            if (free_error == MPI_SUCCESS && layout_error != MPI_SUCCESS)
            {
                free_error = layout_error;
            }
        }
        collective_mpi_check(mpi_comm,
                             free_error,
                             "Failed to free a KS eigenvector MPI file type");

        const int rename_error
            = mpi_rank == 0 && std::rename(temporary_name.c_str(), final_name.c_str()) != 0
                  ? MPI_ERR_IO
                  : MPI_SUCCESS;
        collective_mpi_check(mpi_comm,
                             rename_error,
                             "Failed to publish the completed KS eigenvector MPI file");
    }
    catch (...)
    {
        if (file_opened)
        {
            MPI_File_close(&file);
        }
        for (auto& layout: layouts)
        {
            free_ks_eigenvector_mpi_layout(layout);
        }
        MPI_Barrier(mpi_comm);
        if (mpi_rank == 0)
        {
            std::remove(temporary_name.c_str());
        }
        MPI_Barrier(mpi_comm);
        throw;
    }

    if (mpi_rank == 0)
    {
        std::cout << "KS eigenvector writer: binary v1 MPI-IO, bounded pack buffer, file "
                  << final_name << std::endl;
    }
}
#endif

inline int sum_int_vector(const std::vector<int>& values)
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

template <typename Value>
inline double real_as_double(const Value& value)
{
    return static_cast<double>(value);
}

inline double real_as_double(const std::complex<double>& value)
{
    return value.real();
}

inline bool debug_dump_ewald_split_enabled()
{
    const char* env = std::getenv("ABACUS_DUMP_EWALD_SPLIT_COULOMB");
    if (env == nullptr)
    {
        return false;
    }
    const std::string value(env);
    return !(value.empty() || value == "0" || value == "f" || value == "F"
             || value == "false" || value == "FALSE");
}

inline std::vector<std::vector<int>>
collect_abfs_l_nchi(const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    std::vector<std::vector<int>> abfs_l_nchi;
    abfs_l_nchi.reserve(abfs.size());
    for (const auto& abfs_type : abfs)
    {
        std::vector<int> shell_counts;
        shell_counts.reserve(abfs_type.size());
        for (const auto& abfs_l : abfs_type)
        {
            shell_counts.push_back(static_cast<int>(abfs_l.size()));
        }
        abfs_l_nchi.push_back(std::move(shell_counts));
    }
    return abfs_l_nchi;
}

inline std::vector<std::vector<int>>
collect_wfc_l_nchi(const UnitCell& ucell)
{
    std::vector<std::vector<int>> wfc_l_nchi;
    wfc_l_nchi.reserve(static_cast<std::size_t>(ucell.ntype));
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& atom = ucell.atoms[itype];
        if (atom.nwl < 0 || atom.l_nchi.size() < static_cast<std::size_t>(atom.nwl + 1))
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent AO shell counts.");
        }
        std::vector<int> shell_counts;
        shell_counts.reserve(static_cast<std::size_t>(atom.nwl + 1));
        for (int l = 0; l <= atom.nwl; ++l)
        {
            if (atom.l_nchi[static_cast<std::size_t>(l)] < 0)
            {
                throw std::runtime_error("LibRPA v1 basis output found negative AO shell count.");
            }
            shell_counts.push_back(atom.l_nchi[static_cast<std::size_t>(l)]);
        }
        wfc_l_nchi.push_back(std::move(shell_counts));
    }
    return wfc_l_nchi;
}

inline int basis_size_from_shell_counts(const std::vector<int>& shell_counts)
{
    int basis_size = 0;
    for (std::size_t l = 0; l < shell_counts.size(); ++l)
    {
        const int count = shell_counts[l];
        if (count < 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found negative shell count.");
        }
        const int shell_size = static_cast<int>(2 * l + 1);
        if (count > 0 && shell_size > std::numeric_limits<int>::max() / count)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        const int increment = count * shell_size;
        if (increment > std::numeric_limits<int>::max() - basis_size)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 shell sizes.");
        }
        basis_size += increment;
    }
    return basis_size;
}

inline void write_librpa_split_basis_file(const UnitCell& ucell,
                                          const std::vector<int>& type_sizes,
                                          const std::vector<std::vector<int>>& type_shell_counts,
                                          const std::string& filename)
{
    if (type_sizes.size() != static_cast<std::size_t>(ucell.ntype)
        || type_shell_counts.size() != static_cast<std::size_t>(ucell.ntype))
    {
        throw std::runtime_error("LibRPA v1 basis output found inconsistent atom-type count.");
    }

    int total_basis = 0;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const int type_size = type_sizes[static_cast<std::size_t>(itype)];
        if (type_size <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type basis size.");
        }
        const int shell_size = basis_size_from_shell_counts(type_shell_counts[static_cast<std::size_t>(itype)]);
        if (type_size != shell_size)
        {
            throw std::runtime_error("LibRPA v1 basis output found inconsistent shell layout size.");
        }
        if (ucell.atoms[itype].na < 0 || type_size > std::numeric_limits<int>::max() / ucell.atoms[itype].na)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        const int increment = type_size * ucell.atoms[itype].na;
        if (increment > std::numeric_limits<int>::max() - total_basis)
        {
            throw std::runtime_error("Integer overflow while summing LibRPA v1 basis sizes.");
        }
        total_basis += increment;
    }

    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }
    ofs << std::setw(10) << ucell.ntype
        << std::setw(10) << total_basis
        << "    abacus" << std::endl;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        ofs << std::setw(10) << itype + 1
            << std::setw(10) << type_sizes[static_cast<std::size_t>(itype)]
            << std::endl;
    }
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& shell_counts = type_shell_counts[static_cast<std::size_t>(itype)];
        int nshell = 0;
        for (const int count : shell_counts)
        {
            nshell += count;
        }
        ofs << std::setw(10) << itype + 1
            << std::setw(10) << nshell
            << std::endl;
        for (std::size_t l = 0; l < shell_counts.size(); ++l)
        {
            for (int iradial = 0; iradial < shell_counts[l]; ++iradial)
            {
                ofs << std::setw(10) << l << std::endl;
            }
        }
    }
}

inline void append_unique_abfs_layout_candidates(
    std::vector<std::vector<std::vector<int>>>& candidates_by_type,
    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>& abfs)
{
    if (abfs.empty())
    {
        return;
    }

    const auto shell_counts_by_type = collect_abfs_l_nchi(abfs);
    if (candidates_by_type.empty())
    {
        candidates_by_type.resize(shell_counts_by_type.size());
    }
    else if (candidates_by_type.size() != shell_counts_by_type.size())
    {
        throw std::runtime_error("ABF shell-layout candidates are inconsistent with the atom-type count.");
    }

    for (std::size_t itype = 0; itype < shell_counts_by_type.size(); ++itype)
    {
        const auto& shell_counts = shell_counts_by_type[itype];
        auto& layouts = candidates_by_type[itype];
        if (std::find(layouts.begin(), layouts.end(), shell_counts) == layouts.end())
        {
            layouts.push_back(shell_counts);
        }
    }
}

inline std::vector<std::string> collect_atom_type_labels(const UnitCell& ucell)
{
    std::vector<std::string> labels(static_cast<std::size_t>(ucell.ntype));
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        labels[static_cast<std::size_t>(itype)] = ucell.atoms[itype].label;
    }
    return labels;
}

inline int max_layout_lmax(const std::vector<std::vector<std::vector<int>>>& candidates_by_type)
{
    int lmax = -1;
    for (const auto& type_candidates : candidates_by_type)
    {
        for (const auto& shell_counts : type_candidates)
        {
            lmax = std::max(lmax, static_cast<int>(shell_counts.size()) - 1);
        }
    }
    return lmax;
}

template<typename Tdata>
inline bool has_valid_matrix_shape(const RI::Tensor<Tdata>& tensor)
{
    return tensor.shape.size() == 2 && tensor.shape[0] > 0 && tensor.shape[1] > 0;
}

template<typename Tdata>
inline std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>
collect_local_irreducible_abf_blocks(
    const std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>& period_blocks,
    const std::map<ModuleSymmetry::Tap, std::set<ModuleSymmetry::TC>>& irreducible_sector,
    std::size_t& n_skipped_irreducible_blocks)
{
    std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>> irreducible_blocks;
    n_skipped_irreducible_blocks = 0;
    for (const auto& irap_Rs: irreducible_sector)
    {
        const auto period_iter = period_blocks.find(irap_Rs.first.first);
        for (const auto& irR: irap_Rs.second)
        {
            const RI_2D_Comm::TAC ir_key = {irap_Rs.first.second, irR};
            if (period_iter == period_blocks.end())
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            const auto block_iter = period_iter->second.find(ir_key);
            if (block_iter == period_iter->second.end()
                || !has_valid_matrix_shape(block_iter->second))
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            irreducible_blocks[irap_Rs.first.first][ir_key] = block_iter->second;
        }
    }
    return irreducible_blocks;
}

inline std::size_t sum_skipped_irreducible_blocks(const MPI_Comm& mpi_comm,
                                                  const std::size_t local_count)
{
    unsigned long long global_count = static_cast<unsigned long long>(local_count);
    unsigned long long reduced_count = global_count;
    MPI_Allreduce(&global_count, &reduced_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm);
    return static_cast<std::size_t>(reduced_count);
}

}

template <typename T, typename Tdata>
RPA_LRI<T, Tdata>::~RPA_LRI() = default;

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::postSCF(const UnitCell& ucell,
                                const MPI_Comm& mpi_comm_in,
                                const module_dm::DensityMatrix<T, Tdata>& dm,
                                const elecstate::ElecState* pelec,
                                const K_Vectors& kv,
                                const LCAO_Orbitals& orb,
                                const Parallel_Orbitals& parav,
                                const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("RPA_LRI", "postSCF");
    ModuleBase::timer::start("RPA_LRI", "postSCF");
    ModuleBase::GlobalFunc::MAKE_DIR(outdir);

    if (PARAM.inp.out_librpa_abf_overlap
        && (!PARAM.inp.rpa || PARAM.inp.out_librpa_reader_version != 1
            || this->info.shrink_abfs_pca_thr < 0.0))
    {
        throw std::runtime_error("out_librpa_abf_overlap requires rpa=true, "
                                 "out_librpa_reader_version=1, and shrink ABFs.");
    }

    this->cal_postSCF_exx(dm, mpi_comm_in, ucell, kv, orb, parav);
    if (RpaLriDetail::debug_dump_exx_ao_enabled())
    {
        const std::string file_name_exx
            = PARAM.globalv.global_out_dir + "HexxR" + std::to_string(GlobalV::MY_RANK);
        ModuleIO::write_Hexxs_csr(file_name_exx, ucell, exx_cut_coulomb->Hexxs);
    }
    this->init(mpi_comm_in, kv, orb.cutoffs());
    this->out_bands(pelec);
    this->out_eigen_vector(parav, psi);
    this->out_struc(ucell);
    this->out_bz_sampling();

    std::cout << "rpa_pca_threshold: " << this->info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times: " << this->info.ccp_rmesh_times << std::endl;
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << exx_cut_coulomb->Eexx / 2.0 << std::endl;

    std::cout << "etxc(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha): " << std::fixed << std::setprecision(15)
              << (pelec->f_en.etot - pelec->f_en.etxc + exx_cut_coulomb->Eexx) / 2.0 << std::endl;
    exx_cut_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        cal_large_Cs(ucell, orb, kv);
        cal_abfs_overlap(ucell, orb, kv);
        RpaLriDetail::trim_malloc_cache();
    }
    this->output_ewald_coulomb(ucell, kv, orb);

    ModuleBase::timer::end("RPA_LRI", "postSCF");
}

template <typename T, typename Tdata>
ModuleRI::SternheimerOrbitalSet RPA_LRI<T, Tdata>::take_sternheimer_abfs()
{
    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        return {};
    }
    return std::move(this->abfs);
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::trim_process_heap()
{
    RpaLriDetail::trim_malloc_cache();
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::init(const MPI_Comm& mpi_comm_in, const K_Vectors& kv_in, const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("RPA_LRI", "init");
    ModuleBase::timer::start("RPA_LRI", "init");
    this->mpi_comm = mpi_comm_in;
    this->orb_cutoff_ = orb_cutoff;
    this->lcaos = exx_cut_coulomb->lcaos;
    this->p_kv = &kv_in;
    this->MGT = exx_cut_coulomb->MGT;

    //	this->cv = std::move(exx_lri_rpa.cv);
    //    exx_lri_rpa.cv = exx_lri_rpa.cv;
    ModuleBase::timer::end("RPA_LRI", "init");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_postSCF_exx(const module_dm::DensityMatrix<T, Tdata>& dm,
                                        const MPI_Comm& mpi_comm_in,
                                        const UnitCell& ucell,
                                        const K_Vectors& kv,
                                        const LCAO_Orbitals& orb,
                                        const Parallel_Orbitals& parav)
{
    ModuleBase::TITLE("RPA_LRI", "cal_postSCF_exx");
    ModuleBase::timer::start("RPA_LRI", "cal_postSCF_exx");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv;
    this->orb_cutoff_ = orb.cutoffs();

    Mix_DMk_2D<T> mix_DMk_2D;
    this->use_spacegroup_symmetry_ = (PARAM.inp.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (this->use_spacegroup_symmetry_)
        {mix_DMk_2D.set_nks(kv.get_nkstot_nospin() * (PARAM.inp.nspin == 2 ? 2 : 1));}
    else
        {mix_DMk_2D.set_nks(kv.get_nks());}
        
    // The post-SCF density is only initialized through restart_all().  That
    // path still needs a live mixing engine to initialize per-k-point state;
    // no subsequent mixing step is performed here.
    mix_DMk_2D.set_mixing_plain(0.0);
    if (this->use_spacegroup_symmetry_)
    {
        const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
        const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
        this->symmetry_rotation_.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st, Rs, period, ucell.lat);
        // set Lmax of the rotation matrices to max(l_ao, l_abf), to support rotation under ABF
        this->symmetry_rotation_.set_abfs_Lmax(GlobalC::exx_info.info_ri.abfs_Lmax);
        this->symmetry_rotation_.cal_Ms(kv, ucell, parav, PARAM.inp.nspin);
        mix_DMk_2D.mix(this->symmetry_rotation_.restore_dm(kv, dm.get_dmk_vec(), parav), true);
    }
    else { mix_DMk_2D.mix(dm.get_dmk_vec(), true); }
    
    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>> Ds
        = RI_2D_Comm::split_m2D_ktoR<Tdata>(ucell,
                                            kv,
                                            mix_DMk_2D.get_DMk_out(),
                                            parav,
                                            PARAM.inp.nspin,
                                            this->use_spacegroup_symmetry_);
    
    // reserve exx_ccp_rmesh_times to calculate full Coulomb
    // Note: ccp_type=Hf and hybrid_alpha=1 were previously set on the global Exx_Info
    // and sync_from_global() was called, but this->info (value copy) already has the correct
    // coulomb_param from construction time, so those writes are redundant and removed.
    this->ccp_rmesh_times_ewald = this->info.ccp_rmesh_times;
    // Using rpa_ccp_rmesh_times to calculate cut Coulomb this->Vs_period
    Exx_Info::Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = PARAM.inp.rpa_ccp_rmesh_times;
    if (!exx_cut_coulomb)
        exx_cut_coulomb.reset(new Exx_LRI<double>(local_info));

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);
        this->abfs_shrink = ExxLriDetail::prepare_abfs(
            ucell, orb, this->lcaos, this->info, this->info.shrink_abfs_pca_thr, this->info.files_shrink_abfs);
        const ModuleRI::RpaAbfsPreorthReport preorth_report
            = ModuleRI::finalize_rpa_abfs_from_input(
                this->abfs_shrink, PARAM.inp, PARAM.inp.cal_force);
        if (GlobalV::MY_RANK == 0)
        {
            GlobalV::ofs_running << ModuleRI::format_rpa_abfs_preorth_report(preorth_report);
        }
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_shrink);
    }
    else
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);
        this->abfs = ExxLriDetail::prepare_abfs(
            ucell, orb, this->lcaos, this->info, this->info.pca_threshold, this->info.files_abfs);
        const ModuleRI::RpaAbfsPreorthReport preorth_report
            = ModuleRI::finalize_rpa_abfs_from_input(
                this->abfs, PARAM.inp, PARAM.inp.cal_force);
        if (GlobalV::MY_RANK == 0)
        {
            GlobalV::ofs_running << ModuleRI::format_rpa_abfs_preorth_report(preorth_report);
        }
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, this->abfs);
    }

    if (this->use_spacegroup_symmetry_)
    {
        // Refresh the ABF-side spherical-harmonic rotation matrices after the auxiliary basis
        // is finalized by `init_spencer()`. The earlier `cal_Ms()` call only guaranteed the AO
        // rotation blocks needed for density-matrix restoration.
        this->symmetry_rotation_.set_Cs_rotation(exx_cut_coulomb->get_abfs_nchis());
        this->symmetry_rotation_.cal_Ms(kv, ucell, parav, PARAM.inp.nspin);
    }

    // cal C and V for exx
    this->output_cut_coulomb_cs(ucell, exx_cut_coulomb.get());
    // cal CVCD
    if (this->use_spacegroup_symmetry_ && PARAM.inp.exx_symmetry_realspace)
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav, &this->symmetry_rotation_);
    }
    else
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav);
    }
    // cout<<"postSCF_Eexx: "<<exx_lri_rpa.Eexx<<endl;
    ModuleBase::timer::end("RPA_LRI", "cal_postSCF_exx");
}

// if use shrink, output Coulomb and Cs_data in small abfs
// otherwise, output in normal abfs
template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa)
{
    ModuleBase::TITLE("RPA_LRI", "output_cut_coulomb_cs");
    ModuleBase::timer::start("RPA_LRI", "output_cut_coulomb_cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    std::cout << "Use rpa_ccp_rmesh_times=" << this->info.ccp_rmesh_times << " to calculate cut Coulomb" << std::endl;
    // Shrink_ABFS_ORBITAL cannot exceed this angular momentum of MGT
    exx_lri_rpa->cal_cut_coulomb_cs(Vs_cut_IJR, Cs, ucell, PARAM.inp.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJ
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_cut_IJR, atoms00, atoms01);
    Vs_cut_IJR.clear();
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJ, period);
    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_lri_rpa,
                                  use_shrink ? "basis_aux_shrink_out" : "basis_aux_out",
                                  use_shrink ? "basis_out_shrink" : "basis_out");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "v1_coulomb_cut_iq_", exx_lri_rpa);
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_cut_", exx_lri_rpa);
    }
    Vs_period.clear();
    Vs_period.swap(tmp);

    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_lri_rpa->exx_lri.post_2D.set_tensors_map2(this->Cs_period);

    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
        {
            this->out_Cs_v1(ucell, this->Cs_period, "v1_Cs_shrinked_data_");
        }
        else
        {
            this->out_Cs_v1(ucell, this->Cs_period, "v1_Cs_data_");
        }
    }
    else
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
            this->out_Cs(ucell, this->Cs_period, "Cs_shrinked_data_");
        else
            this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    Cs_period.clear();
    Cs_period.swap(tmp);

    ModuleBase::timer::end("RPA_LRI", "output_cut_coulomb_cs");
}

template <typename T, typename Tdata>
Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<T, Tdata>::select_coulomb_basis_method_(Exx_LRI<double>* exx_lri) const
{
    if (exx_lri == nullptr || exx_lri->exx_objs.empty())
    {
        throw std::invalid_argument("Cannot select Coulomb basis method from an empty Exx_LRI object.");
    }
    return exx_lri->exx_objs.count(Conv_Coulomb_Pot_K::Coulomb_Method::Center2)
        ? Conv_Coulomb_Pot_K::Coulomb_Method::Center2
        : exx_lri->exx_objs.begin()->first;
}

template <typename T, typename Tdata>
std::vector<int> RPA_LRI<T, Tdata>::collect_atom_naux_(const UnitCell& ucell, Exx_LRI<double>* exx_lri) const
{
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> atom_naux(static_cast<std::size_t>(ucell.nat), 0);
    for (int I = 0; I != ucell.nat; ++I)
    {
        atom_naux[static_cast<std::size_t>(I)]
            = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
        if (atom_naux[static_cast<std::size_t>(I)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 output found a non-positive per-atom auxiliary size.");
        }
    }
    return atom_naux;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "output_ewald_coulomb");
    ModuleBase::timer::start("RPA_LRI", "output_ewald_coulomb");

    Exx_Info::Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = this->ccp_rmesh_times_ewald;
    if (!exx_full_coulomb)
        exx_full_coulomb.reset(new Exx_LRI<double>(local_info));

    if (this->info.shrink_abfs_pca_thr >= 0.0)
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs_shrink);
    else
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs);

    const auto write_strict_2d_coulomb_head_sidecar = [&]() {
        if (GlobalC::exx_info.info_ri.ewald_dimension != 2 || GlobalV::MY_RANK != 0)
        {
            return;
        }
        const auto multipoles = Exx_Abfs::Construct_Orbs::get_multipole(exx_full_coulomb->abfs);
        std::vector<std::vector<double>> s_multipoles_by_type(
            static_cast<std::size_t>(ucell.ntype));
        std::vector<int> atoms_per_type(static_cast<std::size_t>(ucell.ntype), 0);
        for (int it = 0; it != ucell.ntype; ++it)
        {
            if (static_cast<std::size_t>(it) < multipoles.size()
                && !multipoles[static_cast<std::size_t>(it)].empty())
            {
                s_multipoles_by_type[static_cast<std::size_t>(it)]
                    = multipoles[static_cast<std::size_t>(it)][0];
            }
            atoms_per_type[static_cast<std::size_t>(it)] = ucell.atoms[it].na;
        }
        const ModuleBase::Vector3<double> a1_bohr = ucell.a1 * ucell.lat0;
        const ModuleBase::Vector3<double> a2_bohr = ucell.a2 * ucell.lat0;
        const auto normalization = RpaLriDetail::strict_2d_coulomb_head_normalization(
            (a1_bohr ^ a2_bohr).norm(), s_multipoles_by_type, atoms_per_type);

        const std::string filename = outdir + "librpa_2d_coulomb_head.dat";
        std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + filename);
        }
        ofs << RpaLriDetail::format_strict_2d_coulomb_head_sidecar(normalization);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to write " + filename);
        }
        std::cout << "Wrote strict 2D Coulomb head normalization to " << filename
                  << ": A_lambda=" << std::setprecision(17)
                  << normalization.raw_head_coefficient
                  << ", sheet_to_raw_scale=" << normalization.sheet_to_raw_scale
                  << std::endl;
    };

    const bool use_direct_2d_coulomb
        = PARAM.inp.out_librpa_2d_coulomb_method == "direct_mixed_fourier";
    const bool use_direct_3d_coulomb
        = PARAM.inp.out_librpa_3d_coulomb_method == "direct_reciprocal";
    if (use_direct_2d_coulomb || use_direct_3d_coulomb)
    {
        if (PARAM.inp.out_librpa_reader_version != 1)
        {
            throw std::invalid_argument(
                "Direct Coulomb output requires out_librpa_reader_version=1.");
        }
        const int required_dimension = use_direct_2d_coulomb ? 2 : 3;
        if (GlobalC::exx_info.info_ri.ewald_dimension != required_dimension)
        {
            throw std::invalid_argument(
                "Direct Coulomb output dimension does not match exx_ewald_dimension.");
        }
        if (RpaLriDetail::debug_dump_ewald_split_enabled()
            || RpaLriDetail::ewald_component_output_enabled())
        {
            throw std::invalid_argument(
                "Direct Coulomb output is incompatible with legacy Ewald split diagnostics.");
        }

        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_full_coulomb.get(),
                                  use_shrink ? "basis_aux_shrink_out" : "basis_aux_out",
                                  use_shrink ? "basis_out_shrink" : "basis_out");
        const auto ewald_object
            = exx_full_coulomb->exx_objs.find(Conv_Coulomb_Pot_K::Coulomb_Method::Ewald);
        if (ewald_object == exx_full_coulomb->exx_objs.end())
        {
            throw std::runtime_error(
                "Direct Coulomb output could not find the Ewald auxiliary object.");
        }
        if (use_direct_2d_coulomb)
        {
            ewald_object->second.evq.output_direct_2d_coulomb(
                ucell,
                PARAM.inp.out_librpa_2d_direct_ecut,
                PARAM.inp.out_librpa_2d_direct_kz_order,
                PARAM.inp.out_librpa_2d_direct_gamma_order);
            write_strict_2d_coulomb_head_sidecar();
        }
        else
        {
            ewald_object->second.evq.output_direct_3d_coulomb(
                ucell, PARAM.inp.out_librpa_3d_direct_ecut);
        }

        exx_full_coulomb.reset();
        RpaLriDetail::trim_malloc_cache();
        ModuleBase::timer::tick("RPA_LRI", "output_ewald_coulomb");
        return;
    }

    // Split Ewald component dumps remain opt-in through the reader-v1 path;
    // keep the legacy split files disabled by default.
    const bool dump_split = false;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    const bool output_ewald_components = RpaLriDetail::ewald_component_output_enabled();
    typename Exx_LRI<double>::EwaldCoulombComponents ewald_components;
    exx_full_coulomb->cal_ewald_coulomb(Vs_full_IJR,
                                        Cs,
                                        ucell,
                                        PARAM.inp.out_ri_cv,
                                        &Vs_short_IJR,
                                        &Vs_long_IJR,
                                        output_ewald_components ? &ewald_components : nullptr);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_ewald, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    const auto gather_periodic = [&](auto& distributed_blocks) {
        auto atom_pair_blocks = RI_2D_Comm::comm_map2_first(mpi_comm, distributed_blocks, atoms00, atoms01);
        distributed_blocks.clear();
        return RI::RI_Tools::cal_period(atom_pair_blocks, period);
    };
    this->Vs_period = gather_periodic(Vs_full_IJR);
    if (this->use_spacegroup_symmetry_)
    {
        // Rebuild the full periodic bare Coulomb from the irreducible real-space sector so that
        // the exported `coulomb_mat_` uses the same ABF rotation convention as the sidecar files.
        this->symmetry_rotation_.set_Cs_rotation(exx_full_coulomb->get_abfs_nchis());

        std::size_t n_skipped_irreducible_blocks_local = 0;
        auto Vs_irreducible = RpaLriDetail::collect_local_irreducible_abf_blocks(
            this->Vs_period,
            this->symmetry_rotation_.get_irreducible_sector(),
            n_skipped_irreducible_blocks_local);
        if (GlobalV::NPROC > 1)
        {
            const std::set<TA> all_atoms_set(atoms.begin(), atoms.end());
            Vs_irreducible = RI_2D_Comm::comm_map2_first(
                this->mpi_comm, Vs_irreducible, all_atoms_set, all_atoms_set);
        }
        const std::size_t n_skipped_irreducible_blocks
            = RpaLriDetail::sum_skipped_irreducible_blocks(
                this->mpi_comm, n_skipped_irreducible_blocks_local);
        if (n_skipped_irreducible_blocks != 0 && GlobalV::MY_RANK == 0)
        {
            std::cout << "Warning: skipped " << n_skipped_irreducible_blocks
                      << " missing irreducible ABF bare-Coulomb blocks during symmetry restoration"
                      << std::endl;
        }
        this->Vs_period
            = this->symmetry_rotation_.restore_HR_abf(ucell.symm, ucell.atoms, ucell.st, Vs_irreducible);
        if (GlobalV::NPROC > 1)
        {
            this->Vs_period = RI_2D_Comm::comm_map2_first(this->mpi_comm, this->Vs_period, atoms00, atoms01);
        }
    }
    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_full_coulomb.get(),
                                  use_shrink ? "basis_aux_shrink_out" : "basis_aux_out",
                                  use_shrink ? "basis_out_shrink" : "basis_out");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "v1_coulomb_full_iq_", exx_full_coulomb.get());
        write_strict_2d_coulomb_head_sidecar();
        if (output_ewald_components)
        {
            auto bare_periodic = gather_periodic(ewald_components.bare_periodic);
            auto gaussian_real = gather_periodic(ewald_components.gaussian_real);
            auto short_range = gather_periodic(ewald_components.short_range);
            auto long_range = gather_periodic(ewald_components.long_range);
            this->out_coulomb_k_v1(
                ucell, bare_periodic, "v1_coulomb_ewald_bare_iq_", exx_full_coulomb.get());
            this->out_coulomb_k_v1(
                ucell, gaussian_real, "v1_coulomb_ewald_gaussian_real_iq_", exx_full_coulomb.get());
            this->out_coulomb_k_v1(
                ucell, short_range, "v1_coulomb_ewald_short_iq_", exx_full_coulomb.get());
            this->out_coulomb_k_v1(
                ucell, long_range, "v1_coulomb_ewald_long_iq_", exx_full_coulomb.get());
        }
    }
    else
    {
        if (output_ewald_components)
        {
            throw std::runtime_error("ABACUS_RPA_EWALD_COMPONENTS requires out_librpa_reader_version=1.");
        }
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_mat_", exx_full_coulomb.get());
    }
    Vs_period.clear();
    Vs_period.swap(tmp);
    Cs.clear();
    Cs.swap(tmp);

    exx_full_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "output_ewald_coulomb");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_large_Cs(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("RPA_LRI", "cal_large_Cs");
    ModuleBase::timer::start("RPA_LRI", "cal_large_Cs");
    if (!exx_cut_coulomb)
        exx_cut_coulomb.reset(new Exx_LRI<double>(this->info));
    this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
    Exx_Abfs::Construct_Orbs::filter_empty_orbs(this->lcaos);
    this->abfs = ExxLriDetail::prepare_abfs(
        ucell, orb, this->lcaos, this->info, this->info.pca_threshold, this->info.files_abfs);
    const ModuleRI::RpaAbfsPreorthReport preorth_report
        = ModuleRI::finalize_rpa_abfs_from_input(
            this->abfs, PARAM.inp, PARAM.inp.cal_force);
    if (GlobalV::MY_RANK == 0)
    {
        GlobalV::ofs_running << ModuleRI::format_rpa_abfs_preorth_report(preorth_report);
    }
    exx_cut_coulomb->init_spencer(this->mpi_comm, ucell, kv, orb, this->abfs);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "exx_cut_coulomb->init");
    this->MGT = exx_cut_coulomb->MGT;
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        atoms[iat] = iat;
    }
    std::map<TA, TatomR> atoms_pos;
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
    const std::array<TatomR, Ndim> latvec = {RI_Util::Vector3_to_array3(ucell.a1),
                                             RI_Util::Vector3_to_array3(ucell.a2),
                                             RI_Util::Vector3_to_array3(ucell.a3)};
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->exx_cut_coulomb->exx_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);
    auto center2_obj_it = this->exx_cut_coulomb->exx_objs.find(Conv_Coulomb_Pot_K::Coulomb_Method::Center2);
    if (center2_obj_it == this->exx_cut_coulomb->exx_objs.end())
    {
        throw std::invalid_argument("RPA_LRI::cal_large_Cs expected a Center2 cut-Coulomb object after init_spencer.");
    }
    center2_obj_it->second.cv.set_orbitals(ucell,
                                           orb,
                                           this->lcaos,
                                           this->abfs,
                                           center2_obj_it->second.abfs_ccp,
                                           this->info.kmesh_times,
                                           this->MGT, // get MGT from exx_cut_coulomb and used in `cal_abfs_overlap`
                                           true);

    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->info.ccp_rmesh_times, ucell, orb_cutoff_);
    std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Vs, 2, false);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Vs start");
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR
        = center2_obj_it->second.cv.cal_Vs(ucell, list_As_Vs.first, list_As_Vs.second[0], {{"writable_Vws", true}});
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Vs end");

    const std::array<Tcell, Ndim> period_Cs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Cs
        = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_Cs, 2, false);
    std::pair<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>,
              std::map<TA, std::map<TAC, std::array<RI::Tensor<Tdata>, 3>>>>
        Cs_dCs = center2_obj_it->second.cv.cal_Cs_dCs(ucell,
                                                      list_As_Cs.first,
                                                      list_As_Cs.second[0],
                                                      {{"cal_dC", false},
                                                       {"writable_Cws", true},
                                                       {"writable_dCws", true},
                                                       {"writable_Vws", false},
                                                       {"writable_dVws", false}});
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_large_Cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    if (PARAM.inp.out_unshrinked_v)
    {
        this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJR, period);
        Vs_cut_IJR.clear();
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Vs_period");
        // MPI: {ia0, {ia1, R}} to {ia0, ia1}
        const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
            = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
        const auto list_A0_pair_R = list_As_Vs_atoms.first;
        const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
        std::set<TA> atoms00;
        std::set<TA> atoms01;
        for (const auto& I: list_A0_pair_R)
        {
            atoms00.insert(I);
        }
        for (const auto& JR: list_A1_pair_R)
        {
            atoms01.insert(JR.first);
        }

        this->Vs_period = RI_2D_Comm::comm_map2_first(this->mpi_comm, this->Vs_period, atoms00, atoms01);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "Vs_period_comm");

        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_unshrinked_cut_", exx_cut_coulomb.get());
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "out_large_Vs");
        this->Vs_period.clear();
        this->Vs_period.swap(tmp);
    }

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs = std::get<0>(Cs_dCs);
    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_cut_coulomb->exx_lri.post_2D.set_tensors_map2(this->Cs_period);
    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        this->out_librpa_basis_v1(ucell, exx_cut_coulomb.get());
        this->out_Cs_v1(ucell, this->Cs_period, "v1_Cs_data_");
    }
    else
    {
        this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "out_large_Cs");
    this->Cs_period.clear();
    this->Cs_period.swap(tmp);
    exx_cut_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    ModuleBase::timer::end("RPA_LRI", "cal_large_Cs");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_abfs_overlap(const UnitCell& ucell, const LCAO_Orbitals& orb, const K_Vectors& kv)
{
    ModuleBase::TITLE("DFT_RPA_interface", "cal_abfs_overlap");
    const auto& abfs_s = this->abfs_shrink;

    // <smaller abfs|smaller abfs>
    Matrix_Orbs11 m_abfs_abfs;
    // <smaller abfs|larger abfs>
    Matrix_Orbs11 m_abfs_abf;

    m_abfs_abf.MGT = this->MGT;
    m_abfs_abf.init(abfs_s, this->abfs, ucell, orb, this->info.kmesh_times);
    m_abfs_abf.init_radial_table();

    m_abfs_abfs.MGT = this->MGT;
    m_abfs_abfs.init(abfs_s, abfs_s, ucell, orb, this->info.kmesh_times);
    m_abfs_abfs.init_radial_table();
    // get Rlist
    const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
    const auto R_period = RI_Util::get_Born_von_Karmen_cells(period);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abfs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abf;

    // index of smaller abfs
    const ModuleBase::Element_Basis_Index::Range range_abfs_s = ModuleBase::Element_Basis_Index::construct_range(abfs_s);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs_s
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs_s);
    // index of larger abfs
    const ModuleBase::Element_Basis_Index::Range range_abfs = ModuleBase::Element_Basis_Index::construct_range(this->abfs);
    const ModuleBase::Element_Basis_Index::IndexLNM index_abfs
        = ModuleBase::Element_Basis_Index::construct_index(range_abfs);

    auto orb_cutoff_ = orb.cutoffs();
    const std::array<Tcell, Ndim> period_Vs = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, std::array<Tcell, Ndim>>>>> list_As_Vs
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);

// Huanjing Gong debug
// std::stringstream ss;
//  ss << "IJR_" << GlobalV::MY_RANK << ".txt";
// std::ofstream ofs;
// ofs.open(ss.str().c_str(), std::ios::out);
// for (size_t iA = 0; iA < list_As_Vs.first.size(); ++iA)
// {
//     const auto& A = list_As_Vs.first[iA];
//     for (const auto& BR: list_As_Vs.second[0])
//     {
//         const auto& B = BR.first;
//         const auto& R = BR.second;
//         ofs << "ABR: " << A << B << "," << R.at(0) << R.at(1) << R.at(2) << std::endl;
//     }
// }
// ofs.close();
#pragma omp parallel
    {
        using LocalMapType = std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>;
        std::map<int, LocalMapType> overlap_abfs_abfs_local;
        std::map<int, LocalMapType> overlap_abfs_abf_local;

#pragma omp for schedule(dynamic) nowait
        for (size_t iA = 0; iA < list_As_Vs.first.size(); ++iA)
        {
            const auto& A = list_As_Vs.first[iA];
            for (const auto& BR: list_As_Vs.second[0])
            {
                const auto& B = BR.first;
                const auto& R = BR.second;

                const size_t TA = ucell.iat2it[A];
                const size_t IA = ucell.iat2ia[A];
                const auto& tauA = ucell.atoms[TA].tau[IA];
                const size_t TB = ucell.iat2it[B];
                const size_t IB = ucell.iat2ia[B];
                const auto& tauB = ucell.atoms[TB].tau[IB];

                const ModuleBase::Vector3<double> tauB_shift
                    = tauB + (RI_Util::array3_to_Vector3(R) * ucell.latvec);
                const ModuleBase::Vector3<double> tau_delta = tauB_shift - tauA;
                static const ModuleBase::Vector3<double> tau0(0.0, 0.0, 0.0);

                auto& local_abfs_abfs = overlap_abfs_abfs_local[A];
                local_abfs_abfs[{B, R}]
                    = m_abfs_abfs.template cal_overlap_matrix<double>(TA,
                                                                      TB,
                                                                      tau0,
                                                                      tau_delta,
                                                                      index_abfs_s,
                                                                      index_abfs_s,
                                                                      Matrix_Orbs11::Matrix_Order::AB);

                auto& local_abfs_abf = overlap_abfs_abf_local[A];
                local_abfs_abf[{B, R}]
                    = m_abfs_abf.template cal_overlap_matrix<double>(TA,
                                                                     TB,
                                                                     tau0,
                                                                     tau_delta,
                                                                     index_abfs_s,
                                                                     index_abfs,
                                                                     Matrix_Orbs11::Matrix_Order::AB);
            }
        }

#pragma omp critical(RPA_LRI_merge)
        {
            for (auto& aPair: overlap_abfs_abfs_local)
            {
                auto& aKey = aPair.first;
                auto& aSubMap = aPair.second;
                for (auto& subPair: aSubMap)
                {
                    auto& key = subPair.first;
                    auto& value = subPair.second;
                    overlap_abfs_abfs[aKey][key] = std::move(value);
                }
            }
            for (auto& aPair: overlap_abfs_abf_local)
            {
                auto& aKey = aPair.first;
                auto& aSubMap = aPair.second;
                for (auto& subPair: aSubMap)
                {
                    auto& key = subPair.first;
                    auto& value = subPair.second;
                    overlap_abfs_abf[aKey][key] = std::move(value);
                }
            }
        }
    }
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    const std::array<Tcell, Ndim> period_Vs_IJ = LRI_CV_Tools::cal_latvec_range<Tcell>(2, ucell, orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abfs_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, overlap_abfs_abfs, atoms00, atoms01);
    overlap_abfs_abfs.clear();
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> overlap_abfs_abf_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, overlap_abfs_abf, atoms00, atoms01);
    overlap_abfs_abf.clear();

    if (PARAM.inp.out_librpa_reader_version == 1)
    {
        out_abfs_overlap_v1(ucell, overlap_abfs_abfs_IJ, overlap_abfs_abf_IJ,
                            "v1_shrink_sinvS_", index_abfs_s, index_abfs);
    }
    else
    {
        out_abfs_overlap(ucell, overlap_abfs_abfs_IJ, overlap_abfs_abf_IJ,
                         "shrink_sinvS_", index_abfs_s, index_abfs);
    }
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                                         std::string filename,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs)
{
    ModuleBase::TITLE("RPA_LRI", "out_abfs_overlap");
    ModuleBase::timer::start("RPA_LRI", "out_abfs_overlap");
    const double threshold = 1e-15;
    const auto format = std::scientific;
    int prec = 15;

    int all_mu_s = 0;
    int all_mu = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    std::vector<int> mu_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        mu_shift[I] = all_mu;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
        all_mu += index_abfs[ucell.iat2it[I]].count_size;
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";

    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;

    // Fourier of ss(R->k), s(R->k)
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_ss;
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_s;
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                RI::Tensor<std::complex<double>> tmp_olp_ss
                    = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                RI::Tensor<std::complex<double>> tmp_olp_s
                    = RI::Global_Func::convert<std::complex<double>>(overlap_abfs_abf[I][{J, R}]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_ss.shape[0], tmp_olp_ss.shape[1]});
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_s.shape[0], tmp_olp_s.shape[1]});
                }
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));

                olp_q_ss[I][{J, q}] = olp_q_ss[I][{J, q}] + tmp_olp_ss * kphase;
                olp_q_s[I][{J, q}] = olp_q_s[I][{J, q}] + tmp_olp_s * kphase;
            }
        }
    }
    // for multi-mpi
    for (int I = 0; I != ucell.nat; I++)
    {
        for (int J = 0; J != ucell.nat; J++)
        {
            for (int ik = 0; ik != nks_tot; ik++)
            {
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs_s[ucell.iat2it[J]].count_size;
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                if (olp_q_s[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs[ucell.iat2it[J]].count_size;
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                for (int ir = 0; ir < olp_q_ss[I][{J, q}].shape[0]; ir++)
                {
                    for (int ic = 0; ic < olp_q_ss[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_ss[I][{J, q}](ir, ic));
                    }
                    for (int ic = 0; ic < olp_q_s[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_s[I][{J, q}](ir, ic));
                    }
                }
            }
        }
    }

    // out_ri_tensor("olp_ss.txt", olp_q_ss, 0.);
    // Inverse of overlap(q)
    inverse_olp(ucell, olp_q_ss, index_abfs_s);
    // out_ri_tensor("olp_ss_inv.txt", olp_q_ss, 0.);
    // out_ri_tensor("olp_s.txt", olp_q_s, 0.);
    for (auto& Ip: overlap_abfs_abf)
    {
        auto I = Ip.first;
        size_t mu_num_s = index_abfs_s[ucell.iat2it[I]].count_size;
        size_t mu_num = index_abfs[ucell.iat2it[I]].count_size;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> sinvS;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                if (sinvS[J].empty())
                {
                    sinvS[J] = RI::Tensor<std::complex<double>>(
                        {overlap_abfs_abfs[I][{J, R}].shape[0], overlap_abfs_abf[I][{J, R}].shape[1]});
                }
            }
            for (const auto& pair: sinvS)
            {
                auto J = pair.first;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                for (int K = 0; K != ucell.nat; K++)
                {
                    sinvS[J] += olp_q_ss.at(I).at({K, q}) * olp_q_s.at(K).at({J, q});
                }
            }
            for (auto& iJU: sinvS)
            {
                auto iJ = iJU.first;
                auto& vq_J = iJU.second;
                size_t nu_num = index_abfs[ucell.iat2it[iJ]].count_size;
                ofs << all_mu_s << "   " << all_mu << "   " << mu_s_shift[I] + 1 << "   " << mu_s_shift[I] + mu_num_s
                    << "  " << mu_shift[iJ] + 1 << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * PARAM.inp.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    // ofs << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].real()
                    //     << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].imag() <<
                    //     std::endl;
                    // if (fabs((*vq_J.data)[i].real()) > threshold || fabs((*vq_J.data)[i].imag()) > threshold)
                    ofs << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].real() << " "
                        << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].imag() << "\n";
                    // else
                    //     ofs << std::showpoint << format << std::setprecision(prec) << 0.0 << " " << std::showpoint
                    //         << format << std::setprecision(prec) << 0.0 << "\n";
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_abfs_overlap");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap_raw_v1(
    const UnitCell& ucell,
    const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    if (!PARAM.inp.rpa || PARAM.inp.out_librpa_reader_version != 1
        || this->info.shrink_abfs_pca_thr < 0.0)
    {
        throw std::runtime_error("raw active-ABF overlap writer requires rpa=true, "
                                 "out_librpa_reader_version=1, and shrink ABFs.");
    }

    const int natom = ucell.nat;
    if (natom <= 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer found no atoms.");
    }
    std::vector<int> atom_naux(static_cast<std::size_t>(natom), 0);
    std::vector<int> atom_shift(static_cast<std::size_t>(natom), 0);
    int naux = 0;
    for (int I = 0; I < natom; ++I)
    {
        const int count = index_abfs_s[ucell.iat2it[I]].count_size;
        if (count <= 0 || naux > std::numeric_limits<int>::max() - count)
        {
            throw std::runtime_error("raw active-ABF overlap writer found an invalid active basis layout.");
        }
        atom_shift[static_cast<std::size_t>(I)] = naux;
        atom_naux[static_cast<std::size_t>(I)] = count;
        naux += count;
    }

    const int nks_tot = PARAM.inp.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    if (nks_tot <= 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer found no q points.");
    }
    const double hermitian_tol = 1e-10;
    const double rank_tol = 1e-12;
    const int natom_metadata = natom;
    int natom_min = 0;
    int natom_max = 0;
    MPI_Allreduce(&natom_metadata, &natom_min, 1, MPI_INT, MPI_MIN, mpi_comm);
    MPI_Allreduce(&natom_metadata, &natom_max, 1, MPI_INT, MPI_MAX, mpi_comm);
    if (natom_min != natom_metadata || natom_max != natom_metadata)
    {
        throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent atom metadata.");
    }
    for (const int count: atom_naux)
    {
        int min_count = 0;
        int max_count = 0;
        MPI_Allreduce(&count, &min_count, 1, MPI_INT, MPI_MIN, mpi_comm);
        MPI_Allreduce(&count, &max_count, 1, MPI_INT, MPI_MAX, mpi_comm);
        if (min_count != count || max_count != count)
        {
            throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent active basis metadata.");
        }
    }

    // comm_map2_first may concentrate all R blocks for one (I,J) on one rank.
    // Check only actual post-communication keys; no complete R coverage is assumed.
    int local_invalid_metadata = 0;
    std::string local_metadata_error;
    std::vector<std::array<int, 5>> local_keys;
    std::vector<int> local_pair_seen(static_cast<std::size_t>(natom) * static_cast<std::size_t>(natom), 0);
    for (const auto& Ip: overlap_abfs_abfs)
    {
        const int I = Ip.first;
        if (I < 0 || I >= natom)
        {
            local_invalid_metadata = 1;
            if (local_metadata_error.empty())
            {
                local_metadata_error = "raw active-ABF overlap writer found an invalid row atom.";
            }
            continue;
        }
        for (const auto& JPp: Ip.second)
        {
            const int J = JPp.first.first;
            const auto R = JPp.first.second;
            if (J < 0 || J >= natom)
            {
                local_invalid_metadata = 1;
                if (local_metadata_error.empty())
                {
                    local_metadata_error = "raw active-ABF overlap writer found invalid block (I,J,R)=("
                                            + std::to_string(I) + "," + std::to_string(J) + ","
                                            + std::to_string(R[0]) + "," + std::to_string(R[1]) + ","
                                            + std::to_string(R[2]) + ").";
                }
                continue;
            }
            const auto& tensor = JPp.second;
            if (tensor.shape.size() != 2
                || tensor.shape[0] != atom_naux[static_cast<std::size_t>(I)]
                || tensor.shape[1] != atom_naux[static_cast<std::size_t>(J)])
            {
                local_invalid_metadata = 1;
                if (local_metadata_error.empty())
                {
                    local_metadata_error = "raw active-ABF overlap writer found inconsistent metadata for block (I,J,R)=("
                                            + std::to_string(I) + "," + std::to_string(J) + ","
                                            + std::to_string(R[0]) + "," + std::to_string(R[1]) + ","
                                            + std::to_string(R[2]) + ").";
                }
                continue;
            }
            local_keys.push_back({I, J, R[0], R[1], R[2]});
            local_pair_seen[static_cast<std::size_t>(I) * static_cast<std::size_t>(natom)
                            + static_cast<std::size_t>(J)] = 1;
        }
    }
    int global_invalid_metadata = 0;
    MPI_Allreduce(&local_invalid_metadata, &global_invalid_metadata, 1, MPI_INT, MPI_MAX, mpi_comm);
    if (global_invalid_metadata != 0)
    {
        throw std::runtime_error(local_metadata_error.empty()
                                     ? "raw active-ABF overlap writer found invalid block metadata on another rank."
                                     : local_metadata_error);
    }
    for (int& seen: local_pair_seen)
    {
        int global_seen = 0;
        MPI_Allreduce(&seen, &global_seen, 1, MPI_INT, MPI_MAX, mpi_comm);
        seen = global_seen;
    }
    for (std::size_t pair_index = 0; pair_index < local_pair_seen.size(); ++pair_index)
    {
        if (local_pair_seen[pair_index] == 0)
        {
            throw std::runtime_error("raw active-ABF overlap writer found a missing atom pair.");
        }
    }

    // Sorting is also a defensive local-map duplicate check; std::map normally
    // makes such a duplicate impossible. The gathered check below rejects both
    // local and cross-rank repeats and reports the same concrete key on all ranks.
    std::sort(local_keys.begin(), local_keys.end());
    const bool local_duplicate = std::adjacent_find(local_keys.begin(), local_keys.end()) != local_keys.end();
    int global_duplicate_flag = 0;
    int mpi_size = 1;
    MPI_Comm_size(mpi_comm, &mpi_size);
    const std::size_t key_width = 5;
    const bool local_count_overflow = local_keys.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / key_width;
    int any_count_overflow = local_count_overflow ? 1 : 0;
    MPI_Allreduce(&any_count_overflow, &global_duplicate_flag, 1, MPI_INT, MPI_MAX, mpi_comm);
    if (global_duplicate_flag != 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer cannot represent MPI key count.");
    }
    const int local_int_count = static_cast<int>(local_keys.size() * key_width);
    std::vector<int> recv_counts(static_cast<std::size_t>(mpi_size), 0);
    MPI_Allgather(&local_int_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, mpi_comm);
    std::vector<int> displacements(static_cast<std::size_t>(mpi_size), 0);
    int total_int_count = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        if (recv_counts[static_cast<std::size_t>(rank)] < 0
            || recv_counts[static_cast<std::size_t>(rank)] > std::numeric_limits<int>::max() - total_int_count)
        {
            throw std::runtime_error("raw active-ABF overlap writer cannot represent gathered MPI key counts.");
        }
        displacements[static_cast<std::size_t>(rank)] = total_int_count;
        total_int_count += recv_counts[static_cast<std::size_t>(rank)];
    }
    std::vector<int> local_packed;
    local_packed.reserve(static_cast<std::size_t>(local_int_count));
    for (const auto& key: local_keys)
    {
        local_packed.insert(local_packed.end(), key.begin(), key.end());
    }
    std::vector<int> gathered(static_cast<std::size_t>(total_int_count));
    MPI_Allgatherv(local_packed.data(), local_int_count, MPI_INT, gathered.data(), recv_counts.data(),
                   displacements.data(), MPI_INT, mpi_comm);
    std::vector<std::array<int, 5>> gathered_keys(static_cast<std::size_t>(total_int_count) / key_width);
    for (std::size_t index = 0; index < gathered_keys.size(); ++index)
    {
        std::copy_n(gathered.begin() + index * key_width, key_width, gathered_keys[index].begin());
    }
    std::sort(gathered_keys.begin(), gathered_keys.end());
    const auto duplicate = std::adjacent_find(gathered_keys.begin(), gathered_keys.end());
    if (local_duplicate || duplicate != gathered_keys.end())
    {
        throw std::runtime_error("raw active-ABF overlap writer found post-communication duplicate contributor for block (I,J,R)=("
                                 + std::to_string((*duplicate)[0]) + "," + std::to_string((*duplicate)[1]) + ","
                                 + std::to_string((*duplicate)[2]) + "," + std::to_string((*duplicate)[3]) + ","
                                 + std::to_string((*duplicate)[4]) + ").");
    }

    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
        const double q_weight = p_kv->wk[ik] / 2.0 * PARAM.inp.nspin;
        if (!std::isfinite(q_weight)
            || !std::isfinite(q[0]) || !std::isfinite(q[1]) || !std::isfinite(q[2]))
        {
            throw std::runtime_error("raw active-ABF overlap writer found non-finite q metadata.");
        }
        for (const double coordinate: q)
        {
            double min_coordinate = 0.0;
            double max_coordinate = 0.0;
            MPI_Allreduce(&coordinate, &min_coordinate, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
            MPI_Allreduce(&coordinate, &max_coordinate, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
            if (max_coordinate - min_coordinate > rank_tol * std::max(1.0, std::abs(coordinate)))
            {
                throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent q coordinates.");
            }
        }
        std::vector<std::complex<double>> overlap(
            static_cast<std::size_t>(naux) * static_cast<std::size_t>(naux), std::complex<double>(0.0, 0.0));
        for (const auto& Ip: overlap_abfs_abfs)
        {
            const int I = Ip.first;
            for (const auto& JPp: Ip.second)
            {
                const int J = JPp.first.first;
                const auto& tensor = JPp.second;
                const auto R = JPp.first.second;
                const double arg = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                    * ModuleBase::TWO_PI;
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                for (int ir = 0; ir < atom_naux[static_cast<std::size_t>(I)]; ++ir)
                {
                    for (int ic = 0; ic < atom_naux[static_cast<std::size_t>(J)]; ++ic)
                    {
                        const std::size_t row = static_cast<std::size_t>(atom_shift[static_cast<std::size_t>(I)] + ir);
                        const std::size_t col = static_cast<std::size_t>(atom_shift[static_cast<std::size_t>(J)] + ic);
                        overlap[row * static_cast<std::size_t>(naux) + col]
                            += static_cast<std::complex<double>>(tensor(ir, ic)) * phase;
                    }
                }
            }
        }
        for (std::complex<double>& value: overlap)
        {
            Parallel_Reduce::reduce_all<std::complex<double>>(value);
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                throw std::runtime_error("raw active-ABF overlap writer found a non-finite overlap value.");
            }
        }
        for (int i = 0; i < naux; ++i)
        {
            for (int j = 0; j < naux; ++j)
            {
                const std::complex<double> difference
                    = overlap[static_cast<std::size_t>(i) * naux + j]
                    - std::conj(overlap[static_cast<std::size_t>(j) * naux + i]);
                const double scale = std::max(
                    1.0,
                    std::max(std::abs(overlap[static_cast<std::size_t>(i) * naux + j]),
                             std::abs(overlap[static_cast<std::size_t>(j) * naux + i])));
                if (std::abs(difference) > hermitian_tol * scale)
                {
                    throw std::runtime_error("raw active-ABF overlap writer found a non-Hermitian S(q).");
                }
            }
        }

        double min_weight = 0.0, max_weight = 0.0;
        MPI_Allreduce(&q_weight, &min_weight, 1, MPI_DOUBLE, MPI_MIN, mpi_comm);
        MPI_Allreduce(&q_weight, &max_weight, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);
        if (max_weight - min_weight > rank_tol * std::max(1.0, std::abs(q_weight)))
        {
            throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent q weight.");
        }

        if (GlobalV::MY_RANK != 0)
        {
            continue;
        }
        const std::string out_name = outdir + "v1_abf_overlap_active_iq_" + std::to_string(ik + 1) + ".dat";
        const std::string tmp_name = out_name + ".tmp";
        std::ofstream ofs(tmp_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + tmp_name);
        }
        const std::int32_t marker = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_MARKER;
        const std::int32_t version = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_VERSION;
        const std::int32_t iq = ik + 1;
        const std::int32_t kind = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_KIND_ACTIVE;
        const std::int32_t naux_i32 = naux;
        const std::int32_t natom_i32 = natom;
        RpaLriDetail::write_scalar(ofs, marker, tmp_name);
        RpaLriDetail::write_scalar(ofs, version, tmp_name);
        RpaLriDetail::write_scalar(ofs, iq, tmp_name);
        RpaLriDetail::write_scalar(ofs, kind, tmp_name);
        RpaLriDetail::write_scalar(ofs, naux_i32, tmp_name);
        RpaLriDetail::write_scalar(ofs, natom_i32, tmp_name);
        RpaLriDetail::write_scalar(ofs, q_weight, tmp_name);
        for (const double coordinate: q)
        {
            RpaLriDetail::write_scalar(ofs, coordinate, tmp_name);
        }
        for (const int count: atom_naux)
        {
            const std::int32_t count_i32 = count;
            RpaLriDetail::write_scalar(ofs, count_i32, tmp_name);
        }
        RpaLriDetail::checked_write(ofs, overlap.data(), overlap.size() * sizeof(std::complex<double>), tmp_name);
        ofs.close();
        if (!ofs.good() || std::rename(tmp_name.c_str(), out_name.c_str()) != 0)
        {
            throw std::runtime_error("Failed to finalize " + out_name);
        }
    }
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap_v1(const UnitCell& ucell,
                                            std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                            std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                                            std::string filename,
                                            const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                                            const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs)
{
    ModuleBase::TITLE("RPA_LRI", "out_abfs_overlap_v1");
    ModuleBase::timer::tick("RPA_LRI", "out_abfs_overlap_v1");

    struct SinvSRecord
    {
        std::int32_t iq = 0;
        std::int32_t nrow_total = 0;
        std::int32_t ncol_total = 0;
        std::int32_t begin_row = 0;
        std::int32_t end_row = 0;
        std::int32_t begin_col = 0;
        std::int32_t end_col = 0;
        double q_weight = 0.0;
        std::int64_t offset = 0;
        std::vector<std::complex<double>> payload;
    };

    int all_mu_s = 0;
    int all_mu = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    std::vector<int> mu_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        mu_shift[I] = all_mu;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
        all_mu += index_abfs[ucell.iat2it[I]].count_size;
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();

    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_ss;
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_s;
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                RI::Tensor<std::complex<double>> tmp_olp_ss
                    = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                RI::Tensor<std::complex<double>> tmp_olp_s
                    = RI::Global_Func::convert<std::complex<double>>(overlap_abfs_abf[I][{J, R}]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_ss.shape[0], tmp_olp_ss.shape[1]});
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_s.shape[0], tmp_olp_s.shape[1]});
                }
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI;
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                olp_q_ss[I][{J, q}] = olp_q_ss[I][{J, q}] + tmp_olp_ss * kphase;
                olp_q_s[I][{J, q}] = olp_q_s[I][{J, q}] + tmp_olp_s * kphase;
            }
        }
    }

    for (int I = 0; I != ucell.nat; I++)
    {
        for (int J = 0; J != ucell.nat; J++)
        {
            for (int ik = 0; ik != nks_tot; ik++)
            {
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs_s[ucell.iat2it[J]].count_size;
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                if (olp_q_s[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs[ucell.iat2it[J]].count_size;
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                for (int ir = 0; ir < olp_q_ss[I][{J, q}].shape[0]; ir++)
                {
                    for (int ic = 0; ic < olp_q_ss[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_ss[I][{J, q}](ir, ic));
                    }
                    for (int ic = 0; ic < olp_q_s[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_s[I][{J, q}](ir, ic));
                    }
                }
            }
        }
    }

    if (PARAM.inp.out_librpa_abf_overlap)
    {
        out_abfs_overlap_raw_v1(ucell, overlap_abfs_abfs, index_abfs_s);
    }

    inverse_olp(ucell, olp_q_ss, index_abfs_s);

    std::vector<SinvSRecord> records;
    for (auto& Ip: overlap_abfs_abf)
    {
        auto I = Ip.first;
        size_t mu_num_s = index_abfs_s[ucell.iat2it[I]].count_size;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> sinvS;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                if (sinvS[J].empty())
                {
                    sinvS[J] = RI::Tensor<std::complex<double>>(
                        {overlap_abfs_abfs[I][{J, R}].shape[0], overlap_abfs_abf[I][{J, R}].shape[1]});
                }
            }
            for (const auto& pair: sinvS)
            {
                auto J = pair.first;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                for (int K = 0; K != ucell.nat; K++)
                {
                    sinvS[J] += olp_q_ss.at(I).at({K, q}) * olp_q_s.at(K).at({J, q});
                }
            }
            for (auto& iJU: sinvS)
            {
                auto iJ = iJU.first;
                auto& vq_J = iJU.second;
                size_t nu_num = index_abfs[ucell.iat2it[iJ]].count_size;
                SinvSRecord record;
                record.iq = static_cast<std::int32_t>(ik + 1);
                record.nrow_total = static_cast<std::int32_t>(all_mu_s);
                record.ncol_total = static_cast<std::int32_t>(all_mu);
                record.begin_row = static_cast<std::int32_t>(mu_s_shift[I] + 1);
                record.end_row = static_cast<std::int32_t>(mu_s_shift[I] + mu_num_s);
                record.begin_col = static_cast<std::int32_t>(mu_shift[iJ] + 1);
                record.end_col = static_cast<std::int32_t>(mu_shift[iJ] + nu_num);
                record.q_weight = p_kv->wk[ik] / 2.0 * PARAM.inp.nspin;
                record.payload.reserve(static_cast<std::size_t>(vq_J.shape[0]) *
                                       static_cast<std::size_t>(vq_J.shape[1]));
                for (int i = 0; i != vq_J.shape[0]; ++i)
                {
                    for (int j = 0; j != vq_J.shape[1]; ++j)
                    {
                        if (i >= static_cast<int>(mu_num_s) || j >= static_cast<int>(nu_num))
                        {
                            throw std::runtime_error("LibRPA v1 shrink_sinvS encountered an inconsistent block shape.");
                        }
                        record.payload.push_back(vq_J(i, j));
                    }
                }
                records.push_back(std::move(record));
            }
        }
    }

    const std::int64_t record_bytes = 7 * static_cast<std::int64_t>(sizeof(std::int32_t))
        + static_cast<std::int64_t>(sizeof(double)) + static_cast<std::int64_t>(sizeof(std::int64_t));
    std::int64_t offset = 2 * static_cast<std::int64_t>(sizeof(std::int32_t))
        + static_cast<std::int64_t>(records.size()) * record_bytes;
    for (auto& record: records)
    {
        record.offset = offset;
        offset += static_cast<std::int64_t>(record.payload.size() * sizeof(std::complex<double>));
    }

    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";
    const std::string out_name = outdir + ss.str();
    std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + out_name);
    }

    const std::int32_t marker = RpaLriDetail::LIBRPA_SHRINK_SINVS_V1_MARKER;
    const std::int32_t nrecords = static_cast<std::int32_t>(records.size());
    RpaLriDetail::write_scalar(ofs, marker, out_name);
    RpaLriDetail::write_scalar(ofs, nrecords, out_name);
    for (const auto& record: records)
    {
        RpaLriDetail::write_scalar(ofs, record.iq, out_name);
        RpaLriDetail::write_scalar(ofs, record.nrow_total, out_name);
        RpaLriDetail::write_scalar(ofs, record.ncol_total, out_name);
        RpaLriDetail::write_scalar(ofs, record.begin_row, out_name);
        RpaLriDetail::write_scalar(ofs, record.end_row, out_name);
        RpaLriDetail::write_scalar(ofs, record.begin_col, out_name);
        RpaLriDetail::write_scalar(ofs, record.end_col, out_name);
        RpaLriDetail::write_scalar(ofs, record.q_weight, out_name);
        RpaLriDetail::write_scalar(ofs, record.offset, out_name);
    }
    for (const auto& record: records)
    {
        RpaLriDetail::checked_write(ofs, record.payload.data(),
                                    record.payload.size() * sizeof(std::complex<double>),
                                    out_name);
    }
    ofs.close();
    ModuleBase::timer::tick("RPA_LRI", "out_abfs_overlap_v1");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::inverse_olp(const UnitCell& ucell,
                                    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                                    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    ModuleBase::TITLE("RPA_LRI", "inverse_olp");
    ModuleBase::timer::start("RPA_LRI", "inverse_olp");
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    size_t all_mu_s = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
    }
    RI::Tensor<std::complex<double>> olp_all = RI::Tensor<std::complex<double>>({all_mu_s, all_mu_s});
    for (int ik = 0; ik < nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                // std::cout << "IJ: " << I << "," << J << std::endl;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;
                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                    {
                        olp_all(mu_s_shift[I] + ir, mu_s_shift[J] + ic) = JPp.second(ir, ic);
                    }
                }
            }
        }
        // for multi-mpi
        // for (int ir = 0; ir < all_mu_s; ir++)
        // {
        //     for (int ic = 0; ic < all_mu_s; ic++)
        //     {
        //         Parallel_Reduce::reduce_all<std::complex<double>>(olp_all(ir, ic));
        //     }
        // }

        // check Hermitian
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                auto delta = std::abs(olp_all(ir, ic) - std::conj(olp_all(ic, ir)));
                if (delta > 1e-10)
                {
                    std::cout << "Warning: olp_all is not Hermitian!" << std::endl;
                    std::cout << "ik,ir,ic: " << ik << "," << ir << "," << ic << std::endl;
                    std::cout << "delta(ir, ic): " << delta << std::endl;
                }
            }
        }
        // out_pure_ri_tensor("olp_all.txt", olp_all, 0.);
        auto olp_inv = LRI_CV_Tools::cal_I(olp_all,
                                           Inverse_Matrix<std::complex<double>>::Method::syev,
                                           this->info.shrink_LU_inv_thr);
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                olp_inv(ic, ir) = std::conj(olp_inv(ir, ic));
            }
        }
        // out_pure_ri_tensor("olp_inv.txt", olp_inv, 0.);
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                auto J = JPp.first.first;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;

                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                        JPp.second(ir, ic) = olp_inv(mu_s_shift[I] + ir, mu_s_shift[J] + ic);
                }
            }
        }
    }
    ModuleBase::timer::end("RPA_LRI", "inverse_olp");
}

// debug function
// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn,
//                                            RI::Tensor<std::complex<double>>& olp,
//                                            const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real()
//                    << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn, RI::Tensor<double>& olp, const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_ri_tensor(const std::string fn,
//                                       std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& olp,
//                                       const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     for (auto& IJq: olp)
//     {
//         int I = IJq.first;
//         for (auto& Jq: IJq.second)
//         {
//             int J = Jq.first.first;
//             auto q = Jq.first.second;
//             auto mat = Jq.second;
//             int nr = mat.shape[0];
//             int nc = mat.shape[1];
//             size_t nnz = nr * nc;
//             fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//             fs << I << " " << J << " " << q.at(0) << " " << q.at(1) << " " << q.at(2) << std::endl;
//             fs << "%" << std::endl;

//             fs << nr << " " << nc << " " << nnz << std::endl;

//             for (int j = 0; j < nc; j++)
//             {
//                 for (int i = 0; i < nr; i++)
//                 {
//                     auto v = mat(i, j);
//                     if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                         fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec)
//                            << v.real() << " " << std::showpoint << format << std::setprecision(prec) << v.imag()
//                            << "\n";
//                 }
//             }
//         }
//     }

//     fs.close();
// }

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = PARAM.inp.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();
    const std::size_t values_per_iw = static_cast<std::size_t>(nbands) * npsin_tmp;

#ifdef __MPI
    const MPI_Comm mpi_comm = parav.comm();
    const int mpi_rank = parav.get_coord_row() * parav.get_dim1() + parav.get_coord_col();
    const int mpi_size = parav.get_dim0() * parav.get_dim1();
    const auto check_mpi = [mpi_comm](const int local_error, const std::string& context) {
        const int local_failed = local_error == MPI_SUCCESS ? 0 : 1;
        int any_failed = 0;
        if (MPI_Allreduce(&local_failed, &any_failed, 1, MPI_INT, MPI_MAX, mpi_comm) != MPI_SUCCESS
            || any_failed != 0)
        {
            throw std::runtime_error(context);
        }
    };
    std::vector<int> output_basis_counts(mpi_size, nbasis / mpi_size);
    for (int ip = 0; ip < nbasis % mpi_size; ++ip)
    {
        ++output_basis_counts[ip];
    }
    std::vector<int> output_basis_offsets(mpi_size + 1, 0);
    std::vector<int> output_basis_owner(nbasis);
    for (int ip = 0; ip < mpi_size; ++ip)
    {
        output_basis_offsets[ip + 1] = output_basis_offsets[ip] + output_basis_counts[ip];
        std::fill(output_basis_owner.begin() + output_basis_offsets[ip],
                  output_basis_owner.begin() + output_basis_offsets[ip + 1], ip);
    }
    const int local_nw = output_basis_counts[mpi_rank];
#else
    const int mpi_rank = 0;
    const int local_nw = nbasis;
#endif
    const std::size_t local_size = static_cast<std::size_t>(local_nw) * values_per_iw;
    std::vector<std::complex<double>> local_wfc(local_size);
#ifdef __MPI
    struct WfcPackIndex
    {
        int spin;
        int band;
        int basis;
    };
    const int local_band_count = parav.ncol_bands;
    const unsigned long long send_size_wide
        = static_cast<unsigned long long>(local_band_count) * psi.get_nbasis() * npsin_tmp;
    const unsigned long long recv_size_wide = static_cast<unsigned long long>(local_nw) * values_per_iw;
    const unsigned long long max_alltoallv_count = std::numeric_limits<int>::max();
    check_mpi(send_size_wide <= max_alltoallv_count && recv_size_wide <= max_alltoallv_count
                  ? MPI_SUCCESS
                  : MPI_ERR_COUNT,
              "RPA eigenvector buffer exceeds MPI_Alltoallv count range.");
    std::vector<int> send_counts(mpi_size, 0);
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (parav.global2local_col(ib) < 0)
            continue;
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
            send_counts[output_basis_owner[parav.local2global_row(ir)]] += npsin_tmp;
    }
    std::vector<int> send_displacements(mpi_size, 0), recv_counts(mpi_size), recv_displacements(mpi_size, 0);
    for (int ip = 1; ip < mpi_size; ++ip)
        send_displacements[ip] = send_displacements[ip - 1] + send_counts[ip - 1];
    check_mpi(MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, mpi_comm),
              "Failed to exchange RPA eigenvector redistribution counts.");
    for (int ip = 1; ip < mpi_size; ++ip)
        recv_displacements[ip] = recv_displacements[ip - 1] + recv_counts[ip - 1];
    const int send_size = static_cast<int>(send_size_wide);
    const int recv_size = static_cast<int>(recv_size_wide);
    check_mpi(recv_displacements.back() + recv_counts.back() == recv_size ? MPI_SUCCESS : MPI_ERR_COUNT,
              "RPA eigenvector Alltoallv redistribution map is inconsistent.");
    std::vector<WfcPackIndex> pack_indices(send_size);
    std::vector<unsigned long long> send_targets(send_size), recv_targets(recv_size);
    std::vector<int> send_positions = send_displacements;
    for (int ib = 0; ib < nbands; ++ib)
    {
        const int ib_local = parav.global2local_col(ib);
        if (ib_local < 0)
            continue;
        for (int ir = 0; ir < psi.get_nbasis(); ++ir)
        {
            const int iw = parav.local2global_row(ir);
            const int destination = output_basis_owner[iw];
            for (int is = 0; is < npsin_tmp; ++is)
            {
                const int position = send_positions[destination]++;
                pack_indices[position] = {is, ib_local, ir};
                send_targets[position] = static_cast<unsigned long long>(iw - output_basis_offsets[destination])
                    * values_per_iw + ib * npsin_tmp + is;
            }
        }
    }
    std::complex<double> dummy(0.0, 0.0);
    unsigned long long dummy_target = 0;
    check_mpi(MPI_Alltoallv(send_size > 0 ? send_targets.data() : &dummy_target,
                            send_counts.data(), send_displacements.data(), MPI_UNSIGNED_LONG_LONG,
                            recv_size > 0 ? recv_targets.data() : &dummy_target,
                            recv_counts.data(), recv_displacements.data(), MPI_UNSIGNED_LONG_LONG, mpi_comm),
              "Failed to exchange RPA eigenvector destination indices.");
    std::vector<std::complex<double>> send_wfc(send_size), recv_wfc(recv_size);
#endif

    if (PARAM.inp.out_librpa_reader_version == 1)
    {
#ifdef __MPI
        ModuleBase::timer::tick("RPA_LRI", "out_eigen_vector_v1_mpi_io");
        try
        {
            const MPI_Comm io_comm = parav.comm();
            if (io_comm == MPI_COMM_NULL)
            {
                throw std::runtime_error("KS eigenvector MPI-IO writer has no wavefunction communicator.");
            }
            int communicator_relation = MPI_UNEQUAL;
            RpaLriDetail::collective_mpi_check(
                io_comm,
                MPI_Comm_compare(io_comm, this->mpi_comm, &communicator_relation),
                "Failed to compare KS eigenvector MPI communicators");
            RpaLriDetail::collective_require(
                io_comm,
                communicator_relation == MPI_IDENT || communicator_relation == MPI_CONGRUENT,
                "KS eigenvector wavefunction and RPA communicators are inconsistent");
            RpaLriDetail::write_ks_eigenvector_v1_mpi(io_comm,
                                                      parav,
                                                      psi,
                                                      nks_tot,
                                                      npsin_tmp,
                                                      PARAM.inp.nspin,
                                                      PARAM.inp.nbands,
                                                      PARAM.globalv.nlocal,
                                                      outdir + "KS_eigenvector_0.dat");
        }
        catch (...)
        {
            ModuleBase::timer::tick("RPA_LRI", "out_eigen_vector_v1_mpi_io");
            throw;
        }
        ModuleBase::timer::tick("RPA_LRI", "out_eigen_vector_v1_mpi_io");
#else
        struct KSEigenRecord
        {
            std::int32_t ik = 0;
            std::int64_t payload_offset = 0;
            std::vector<std::complex<double>> payload;
        };

        std::vector<KSEigenRecord> records;
        records.reserve(static_cast<std::size_t>(nks_tot));

        for (int ik = 0; ik < nks_tot; ik++)
        {
            std::vector<ModuleBase::ComplexMatrix> is_wfc_ib_iw(npsin_tmp);
            for (int is = 0; is < npsin_tmp; is++)
            {
                is_wfc_ib_iw[is].create(PARAM.inp.nbands, PARAM.globalv.nlocal);
                for (int ib_global = 0; ib_global < PARAM.inp.nbands; ++ib_global)
                {
                    std::vector<std::complex<double>> wfc_iks(PARAM.globalv.nlocal, zero);

                    const int ib_local = parav.global2local_col(ib_global);

                    if (ib_local >= 0)
                    {
                        for (int ir = 0; ir < psi.get_nbasis(); ir++)
                        {
                            wfc_iks[parav.local2global_row(ir)] = psi(ik + nks_tot * is, ib_local, ir);
                        }
                    }

                    for (int iw = 0; iw < PARAM.globalv.nlocal; iw++)
                    {
                        is_wfc_ib_iw[is](ib_global, iw) = wfc_iks[iw];
                    }
                }
            }

            if (GlobalV::MY_RANK == 0)
            {
                KSEigenRecord record;
                record.ik = static_cast<std::int32_t>(ik + 1);
                record.payload.reserve(static_cast<std::size_t>(npsin_tmp)
                                       * static_cast<std::size_t>(PARAM.inp.nbands)
                                       * static_cast<std::size_t>(PARAM.globalv.nlocal));
                if (PARAM.inp.nspin == 4)
                {
                    if (PARAM.globalv.nlocal % 2 != 0)
                    {
                        throw std::runtime_error("SOC KS eigenvector output expects an even basis size.");
                    }
                    const int nlocal_ao = PARAM.globalv.nlocal / 2;
                    for (int isoc = 0; isoc < 2; ++isoc)
                    {
                        for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
                        {
                            for (int iw = 0; iw < nlocal_ao; ++iw)
                            {
                                record.payload.push_back(is_wfc_ib_iw[0](ib, iw * 2 + isoc));
                            }
                        }
                    }
                }
                else
                {
                    for (int is = 0; is < npsin_tmp; ++is)
                    {
                        for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
                        {
                            for (int iw = 0; iw < PARAM.globalv.nlocal; ++iw)
                            {
                                record.payload.push_back(is_wfc_ib_iw[is](ib, iw));
                            }
                        }
                    }
                }
                records.push_back(std::move(record));
            }
        }

        if (GlobalV::MY_RANK == 0)
        {
            const std::string out_name = outdir + "KS_eigenvector_0.dat";
            const std::int64_t record_bytes = static_cast<std::int64_t>(sizeof(std::int32_t))
                + static_cast<std::int64_t>(sizeof(std::int64_t));
            std::int64_t offset = 6 * static_cast<std::int64_t>(sizeof(std::int32_t))
                + static_cast<std::int64_t>(records.size()) * record_bytes;
            for (auto& record: records)
            {
                record.payload_offset = offset;
                offset += static_cast<std::int64_t>(record.payload.size() * sizeof(std::complex<double>));
            }

            std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
            if (!ofs.good())
            {
                throw std::runtime_error("Failed to open " + out_name);
            }

            const std::int32_t marker = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_MARKER;
            const std::int32_t kind = RpaLriDetail::LIBRPA_KS_EIGENVECTOR_V1_KIND_COMPLEX_DOUBLE;
            const std::int32_t nkpoints_local = RpaLriDetail::checked_i32_from_size(records.size(),
                                                                                    "KS eigenvector k-point count");
            const std::int32_t nspins = RpaLriDetail::checked_i32_from_int(npsin_tmp, "KS eigenvector spin count");
            const std::int32_t nstates = RpaLriDetail::checked_i32_from_int(PARAM.inp.nbands,
                                                                            "KS eigenvector state count");
            const std::int32_t nbasis_wfc = RpaLriDetail::checked_i32_from_int(PARAM.globalv.nlocal,
                                                                               "KS eigenvector basis count");
            RpaLriDetail::write_scalar(ofs, marker, out_name);
            RpaLriDetail::write_scalar(ofs, kind, out_name);
            RpaLriDetail::write_scalar(ofs, nkpoints_local, out_name);
            RpaLriDetail::write_scalar(ofs, nspins, out_name);
            RpaLriDetail::write_scalar(ofs, nstates, out_name);
            RpaLriDetail::write_scalar(ofs, nbasis_wfc, out_name);
            for (const auto& record: records)
            {
                RpaLriDetail::write_scalar(ofs, record.ik, out_name);
                RpaLriDetail::write_scalar(ofs, record.payload_offset, out_name);
            }
            for (const auto& record: records)
            {
                RpaLriDetail::checked_write(ofs,
                                            record.payload.data(),
                                            record.payload.size() * sizeof(std::complex<double>),
                                            out_name);
            }
            ofs.close();
        }
#endif
        return;
    }

    std::string output_buffer;
    constexpr std::size_t output_line_bytes = 61;
    output_buffer.reserve(local_size * output_line_bytes + 32);
    const std::string filename = outdir + "KS_eigenvector.txt";
#ifdef __MPI
    MPI_File file = MPI_FILE_NULL;
    unsigned long long total_file_bytes = 0;
    const unsigned long long max_mpi_offset = static_cast<unsigned long long>(std::numeric_limits<MPI_Offset>::max());
    const char dummy_buffer = '\0';
    check_mpi(MPI_File_open(mpi_comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY,
                            MPI_INFO_NULL, &file), "Failed to open " + filename + ".");
    check_mpi(MPI_File_set_size(file, 0), "Failed to truncate " + filename + ".");
#else
    std::ofstream ofs(filename.c_str(), std::ios::out);
#endif

    for (int ik = 0; ik < nks_tot; ik++)
    {
#ifdef __MPI
        // 4.1 Pack and redistribute all bands and spins for this k point.
        for (int index = 0; index < send_size; ++index)
        {
            const WfcPackIndex& source = pack_indices[index];
            send_wfc[index] = psi(ik + nks_tot * source.spin, source.band, source.basis);
        }
        check_mpi(MPI_Alltoallv(send_size > 0 ? send_wfc.data() : &dummy,
                                send_counts.data(),
                                send_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                recv_size > 0 ? recv_wfc.data() : &dummy,
                                recv_counts.data(),
                                recv_displacements.data(),
                                MPI_DOUBLE_COMPLEX,
                                mpi_comm),
                  "Failed to redistribute RPA eigenvectors with MPI_Alltoallv.");
        for (int index = 0; index < recv_size; ++index)
        {
            local_wfc[static_cast<std::size_t>(recv_targets[index])] = recv_wfc[index];
        }
#else
        for (int is = 0; is < npsin_tmp; is++)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int iw = 0; iw < psi.get_nbasis(); ++iw)
                {
                    local_wfc[iw * values_per_iw + ib * npsin_tmp + is]
                        = psi(ik + nks_tot * is, ib, iw);
                }
            }
        }
#endif
        // 4.2 Format this rank's contiguous basis interval in basis-band-spin order.
        output_buffer.clear();
        if (mpi_rank == 0)
        {
            output_buffer += std::to_string(ik + 1);
            output_buffer += '\n';
        }
        char line_buffer[128];
        for (int iw = 0; iw < local_nw; ++iw)
        {
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int is = 0; is < npsin_tmp; is++)
                {
                    const std::complex<double>& value
                        = local_wfc[static_cast<std::size_t>(iw) * values_per_iw + ib * npsin_tmp + is];
                    const int line_length = std::snprintf(line_buffer,
                                                          sizeof(line_buffer),
                                                          "%30.15f%30.15f\n",
                                                          value.real(),
                                                          value.imag());
                    if (line_length != static_cast<int>(output_line_bytes))
                    {
                        throw std::runtime_error("Failed to format an RPA eigenvector value.");
                    }
                    output_buffer.append(line_buffer, static_cast<std::size_t>(line_length));
                }
            }
        }

#ifdef __MPI
        // 4.3 Fixed-width value lines make each rank's file offset deterministic;
        // only the k-point header length varies.
        const unsigned long long local_bytes
            = static_cast<unsigned long long>(output_buffer.size());
        const unsigned long long header_bytes = std::to_string(ik + 1).size() + 1;
        const unsigned long long bytes_per_basis = values_per_iw * output_line_bytes;
        if (bytes_per_basis > (max_mpi_offset - header_bytes) / std::max(nbasis, 1)
            || header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis
                   > max_mpi_offset - total_file_bytes)
        {
            throw std::runtime_error("RPA eigenvector file layout exceeds MPI_Offset range.");
        }
        const unsigned long long kpoint_bytes
            = header_bytes + static_cast<unsigned long long>(nbasis) * bytes_per_basis;
        const unsigned long long rank_offset_bytes
            = mpi_rank == 0
                  ? 0
                  : header_bytes
                        + static_cast<unsigned long long>(output_basis_offsets[mpi_rank])
                              * bytes_per_basis;
        const unsigned long long file_offset_bytes = total_file_bytes + rank_offset_bytes;
        // Rank 0 has the largest basis slice and additionally owns the k-point header.
        const unsigned long long max_local_bytes
            = header_bytes
              + static_cast<unsigned long long>(output_basis_counts[0]) * bytes_per_basis;
        const unsigned long long max_write_count
            = static_cast<unsigned long long>(std::numeric_limits<int>::max());
        // MPI_File_write_at_all is collective, so every rank executes the same
        // number of chunks; ranks without data in a chunk pass a zero count.
        for (unsigned long long buffer_offset = 0;
             buffer_offset < max_local_bytes;
             buffer_offset += max_write_count)
        {
            const unsigned long long remaining
                = local_bytes > buffer_offset ? local_bytes - buffer_offset : 0;
            const int write_count
                = static_cast<int>(std::min(remaining, max_write_count));
            const char* write_buffer
                = write_count > 0
                      ? output_buffer.data() + static_cast<std::size_t>(buffer_offset)
                      : &dummy_buffer;
            const unsigned long long write_offset
                = write_count > 0 ? file_offset_bytes + buffer_offset : file_offset_bytes;
            check_mpi(MPI_File_write_at_all(file,
                                            static_cast<MPI_Offset>(write_offset),
                                            write_buffer,
                                            write_count,
                                            MPI_CHAR,
                                            MPI_STATUS_IGNORE),
                      "Failed to write " + filename + ".");
        }
        total_file_bytes += kpoint_bytes;
#else
        ofs << output_buffer;
#endif
    } // ik
#ifdef __MPI
    check_mpi(MPI_File_sync(file), "Failed to sync " + filename + ".");
    check_mpi(MPI_File_close(&file), "Failed to close " + filename + ".");
#endif
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_struc(const UnitCell& ucell)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    const auto unit_scales = RpaLriDetail::librpa_stru_unit_scales(ucell.lat0);
    const ModuleBase::Matrix3 lat = ucell.latvec * unit_scales.real_space_bohr;
    const ModuleBase::Matrix3 G_RPA = ucell.G * unit_scales.reciprocal_space_bohr_inv;
    std::ofstream ofs;
    ofs.open(outdir + "stru_out.txt", std::ios::out);
    const auto write_scientific_triplet = [&ofs](const double x, const double y, const double z) {
        ofs << std::setw(24) << std::scientific << std::setprecision(15) << x
            << std::setw(24) << std::scientific << std::setprecision(15) << y
            << std::setw(24) << std::scientific << std::setprecision(15) << z << std::endl;
    };
    write_scientific_triplet(lat.e11, lat.e12, lat.e13);
    write_scientific_triplet(lat.e21, lat.e22, lat.e23);
    write_scientific_triplet(lat.e31, lat.e32, lat.e33);

    write_scientific_triplet(G_RPA.e11, G_RPA.e12, G_RPA.e13);
    write_scientific_triplet(G_RPA.e21, G_RPA.e22, G_RPA.e23);
    write_scientific_triplet(G_RPA.e31, G_RPA.e32, G_RPA.e33);

    ofs << ucell.nat << std::endl;
    for (int it = 0; it < ucell.ntype; it++)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const auto position_bohr = ucell.atoms[it].tau[ia] * unit_scales.real_space_bohr;
            ofs << std::setw(24) << std::scientific << std::setprecision(15) << position_bohr.x
                << std::setw(24) << std::scientific << std::setprecision(15) << position_bohr.y
                << std::setw(24) << std::scientific << std::setprecision(15) << position_bohr.z
                << std::setw(15) << (it + 1) << std::endl;
        }
    }

    if (ModuleSymmetry::Symmetry::symm_flag == 1 && ucell.symm.nrotk > 0)
    {
        const auto& symm = ucell.symm;
        // Export the spatial parts of the symmetry operations in one common row block.
        const int n_anti = symm.magnetic_nspin4 ? symm.nrotk_anti : 0;
        std::vector<RpaLriDetail::LibRpaSymmetryOperation> unitary(symm.nrotk);
        for (int isym = 0; isym < symm.nrotk; ++isym)
        {
            unitary[isym] = RpaLriDetail::make_librpa_symmetry_operation(symm.gmatrix[isym], symm.gtrans[isym]);
        }
        std::vector<RpaLriDetail::LibRpaSymmetryOperation> antiunitary(n_anti);
        for (int isym = 0; isym < n_anti; ++isym)
        {
            antiunitary[isym]
                = RpaLriDetail::make_librpa_symmetry_operation(symm.gmatrix_anti[isym], symm.gtrans_anti[isym]);
        }
        RpaLriDetail::write_librpa_symmetry_rows(ofs, unitary, antiunitary);

        if (PARAM.inp.nspin == 4)
        {
            // LibRPA consumes the explicit spin-space table for SOC. ABACUS stores
            // symmetry rotations in fractional row-vector form, so convert each
            // operation to Cartesian coordinates before constructing its SU(2) part.
            const auto spin_operation = [&ucell](const ModuleBase::Matrix3& rotation,
                                                 const int antiunitary) {
                const ModuleBase::Matrix3 cartesian = ucell.latvec.Inverse() * rotation * ucell.latvec;
                const auto spin_u = ModuleSymmetry::SpinRotation::so3_to_su2(cartesian);
                RpaLriDetail::LibRpaSpinSymmetryOperation operation;
                operation.antiunitary = antiunitary;
                for (int i = 0; i != 4; ++i)
                {
                    operation.spin_u[static_cast<std::size_t>(2 * i)] = spin_u[i].real();
                    operation.spin_u[static_cast<std::size_t>(2 * i + 1)] = spin_u[i].imag();
                }
                return operation;
            };
            std::vector<RpaLriDetail::LibRpaSpinSymmetryOperation> spin_operations;
            spin_operations.reserve(static_cast<std::size_t>(symm.nrotk + n_anti));
            for (int isym = 0; isym != symm.nrotk; ++isym)
            {
                spin_operations.push_back(spin_operation(symm.gmatrix[isym], 0));
            }
            for (int isym = 0; isym != n_anti; ++isym)
            {
                spin_operations.push_back(spin_operation(symm.gmatrix_anti[isym], 1));
            }
            const int grey_group = symm.magnetic_nspin4 ? 0 : 1;
            RpaLriDetail::write_librpa_spin_symmetry(ofs, grey_group, 1, spin_operations);
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bz_sampling()
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    ModuleBase::TITLE("DFT_RPA_interface", "out_bz_sampling");
    const int nks_tot = PARAM.inp.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const int n_coulomb_irreducible = RpaLriDetail::librpa_stored_coulomb_q_count(nks_tot);

    const std::string filename = outdir + "bz_sampling_out";
    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }
    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::endl;
    ofs << nks_tot << std::setw(8) << n_coulomb_irreducible << std::endl;
    double weight_sum = 0.0;
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        weight_sum += p_kv->wk[ik];
    }
    if (weight_sum <= 0.0)
    {
        throw std::runtime_error("Cannot write " + filename + " with non-positive total k-point weight.");
    }
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const auto stored_q_index = RpaLriDetail::librpa_stored_q_index(ik, nks_tot);
        ofs << std::setw(8) << ik + 1
            << std::setw(24) << std::scientific << std::setprecision(15)
            << (p_kv->wk[ik] / weight_sum)
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].x
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].y
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].z
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_c[ik].x
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_c[ik].y
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_c[ik].z
            << std::setw(8) << stored_q_index.coulomb_irreducible_index
            << std::setw(8) << stored_q_index.representative_scf_index
            << std::endl;
    }
    ofs.close();
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bands(const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    std::ofstream ofs;
    ofs.open(outdir + "band_out.txt", std::ios::out);
    // Set precision before the Fermi energy and first occupation are written.
    // The occupations include k weights; rounding the first row changes the
    // normalized occupation when LibRPA restores an irreducible k mesh.
    ofs << std::fixed << std::setprecision(15);
    ofs << nks_tot << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << PARAM.inp.nbands << std::endl;
    ofs << PARAM.globalv.nlocal << std::endl;
    ofs << (pelec->eferm.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != PARAM.inp.nbands; ib++)
            {
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << pelec->wg(ik + is * nks_tot, ib) * nks_tot
                    << std::setw(25) << pelec->ekb(ik + is * nks_tot, ib) / 2.0
                    << std::setw(25)
                    << pelec->ekb(ik + is * nks_tot, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs");
    ModuleBase::timer::start("RPA_LRI", "out_Cs");

    std::stringstream ss;
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";
    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);
    ofs << ucell.nat << "    " << 0 << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Cs_in)
    {
        size_t I = Ip.first;
        size_t i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            size_t J = JPp.first.first;
            auto R = JPp.first.second;
            auto& tmp_Cs = JPp.second;
            size_t j_num = ucell.atoms[ucell.iat2it[J]].nw;

            ofs << I + 1 << "   " << J + 1 << "   " << R[0] << "   " << R[1] << "   " << R[2] << "   " << i_num
                << std::endl;
            ofs << j_num << "   " << tmp_Cs.shape[0] << std::endl;
            for (int i = 0; i != i_num; i++)
            {
                for (int j = 0; j != j_num; j++)
                {
                    for (int mu = 0; mu != tmp_Cs.shape[0]; mu++)
                    {
                        ofs << std::setw(30) << tmp_Cs(mu, i, j) << std::endl;
                    }
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_Cs");
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs_v1(const UnitCell& ucell,
                                  std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in,
                                  std::string filename)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs_v1");
    ModuleBase::timer::tick("RPA_LRI", "out_Cs_v1");

    struct CsRecord
    {
        int ia1 = 0;
        int ia2 = 0;
        int R[3] = {0, 0, 0};
        double max_abs = 0.0;
        std::int64_t offset = 0;
        RI::Tensor<Tdata>* tensor = nullptr;
        int nw1 = 0;
        int nw2 = 0;
        int naux = 0;
    };

    std::vector<CsRecord> records;
    records.reserve(Cs_in.size());
    for (auto& Ip: Cs_in)
    {
        const int I = static_cast<int>(Ip.first);
        const int i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            const int J = static_cast<int>(JPp.first.first);
            auto& tmp_Cs = JPp.second;
            if (tmp_Cs.shape.size() != 3 || tmp_Cs.shape[0] <= 0 || tmp_Cs.shape[1] <= 0 || tmp_Cs.shape[2] <= 0)
            {
                continue;
            }
            const int j_num = ucell.atoms[ucell.iat2it[J]].nw;
            if (static_cast<int>(tmp_Cs.shape[1]) != i_num || static_cast<int>(tmp_Cs.shape[2]) != j_num)
            {
                throw std::runtime_error("LibRPA v1 Cs output encountered an inconsistent tensor shape.");
            }
            CsRecord record;
            record.ia1 = I + 1;
            record.ia2 = J + 1;
            record.R[0] = JPp.first.second[0];
            record.R[1] = JPp.first.second[1];
            record.R[2] = JPp.first.second[2];
            record.tensor = &tmp_Cs;
            record.nw1 = i_num;
            record.nw2 = j_num;
            record.naux = static_cast<int>(tmp_Cs.shape[0]);
            for (int i = 0; i != record.nw1; ++i)
            {
                for (int j = 0; j != record.nw2; ++j)
                {
                    for (int mu = 0; mu != record.naux; ++mu)
                    {
                        record.max_abs = std::max(record.max_abs, std::abs(RpaLriDetail::real_as_double(tmp_Cs(mu, i, j))));
                    }
                }
            }
            records.push_back(record);
        }
    }

    const std::int64_t nblocks = static_cast<std::int64_t>(records.size());
    const std::int64_t record_bytes = static_cast<std::int64_t>(5 * sizeof(std::int32_t)
        + sizeof(double) + sizeof(std::int64_t));
    std::int64_t offset = static_cast<std::int64_t>(3 * sizeof(std::int32_t) + 2 * sizeof(std::int64_t))
        + nblocks * record_bytes;
    for (auto& record: records)
    {
        record.offset = offset;
        const unsigned long long nw_product = RpaLriDetail::checked_mul_u64(
            static_cast<unsigned long long>(record.nw1),
            static_cast<unsigned long long>(record.nw2),
            "LibRPA v1 Cs block size");
        const unsigned long long values = RpaLriDetail::checked_mul_u64(
            nw_product,
            static_cast<unsigned long long>(record.naux),
            "LibRPA v1 Cs block size");
        const unsigned long long bytes = RpaLriDetail::checked_mul_u64(
            values,
            static_cast<unsigned long long>(sizeof(double)),
            "LibRPA v1 Cs block size");
        offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Cs block size");
    }

    std::stringstream ss;
    ss << filename << GlobalV::MY_RANK << ".txt";
    const std::string out_name = outdir + ss.str();
    std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + out_name);
    }

    const std::int32_t marker = RpaLriDetail::LIBRPA_LRICOEF_V1_MARKER;
    const std::int32_t natom = static_cast<std::int32_t>(ucell.nat);
    const std::int32_t ncell = 0;
    RpaLriDetail::write_scalar(ofs, marker, out_name);
    RpaLriDetail::write_scalar(ofs, natom, out_name);
    RpaLriDetail::write_scalar(ofs, ncell, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    for (const auto& record: records)
    {
        const std::int32_t ia1 = record.ia1;
        const std::int32_t ia2 = record.ia2;
        const std::int32_t r0 = record.R[0];
        const std::int32_t r1 = record.R[1];
        const std::int32_t r2 = record.R[2];
        RpaLriDetail::write_scalar(ofs, ia1, out_name);
        RpaLriDetail::write_scalar(ofs, ia2, out_name);
        RpaLriDetail::write_scalar(ofs, r0, out_name);
        RpaLriDetail::write_scalar(ofs, r1, out_name);
        RpaLriDetail::write_scalar(ofs, r2, out_name);
        RpaLriDetail::write_scalar(ofs, record.max_abs, out_name);
        RpaLriDetail::write_scalar(ofs, record.offset, out_name);
    }
    for (const auto& record: records)
    {
        for (int i = 0; i != record.nw1; ++i)
        {
            for (int j = 0; j != record.nw2; ++j)
            {
                for (int mu = 0; mu != record.naux; ++mu)
                {
                    const double value = RpaLriDetail::real_as_double((*record.tensor)(mu, i, j));
                    RpaLriDetail::write_scalar(ofs, value, out_name);
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::tick("RPA_LRI", "out_Cs_v1");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k(const UnitCell& ucell,
                                      std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                      std::string filename,
                                      Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k");
    ModuleBase::timer::start("RPA_LRI", "out_coulomb_k");

    int all_mu = 0;
    std::vector<int> mu_shift(ucell.nat);
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
    }
    const int nks_tot = PARAM.inp.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (GlobalV::MY_RANK + 1) << ".txt";

    std::ofstream ofs;
    ofs.open(outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Vs)
    {
        auto I = Ip.first;
        size_t mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;

                auto R = JPp.first.second;
                if (J < I)
                {
                    continue;
                }
                if (!RpaLriDetail::has_valid_matrix_shape(JPp.second))
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                auto iJ = vq_Jp.first;
                auto& vq_J = vq_Jp.second;
                size_t nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[iJ]);
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[iJ] + 1
                    << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * PARAM.inp.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    ofs << std::setw(25) << (*vq_J.data)[i].real()
                        << std::setw(25) << (*vq_J.data)[i].imag() << std::endl;
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_coulomb_k");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_velocity(const UnitCell& ucell,
                                     const Grid_Driver& gd,
                                     const TwoCenterBundle& two_center_bundle,
                                     const Parallel_Orbitals& parav,
                                     const psi::Psi<T>& psi,
                                     const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_velocity");
    ModuleBase::timer::start("RPA_LRI", "out_velocity");

    Parallel_2D parac;
    LR_Util::setup_2d_division(parac, parav.get_block_size(), PARAM.globalv.nlocal, PARAM.inp.nbands
#ifdef __MPI
                               , parav.blacs_ctxt
#endif
    );

    const int nk = PARAM.inp.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = PARAM.inp.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();

    std::vector<int> nocc(2, nbands);
    std::vector<int> nvirt(2, 0);
    const std::vector<std::complex<double>> velocity_mo
        = LR_Util::cal_velocity_mo(ucell,
                                   gd,
                                   two_center_bundle,
                                   parav,
                                   parac,
                                   *this->p_kv,
                                   psi,
                                   nk,
                                   nspin_tmp,
                                   PARAM.globalv.nlocal,
                                   nocc,
                                   nvirt);
    if (GlobalV::MY_RANK == 0)
    {
        LR_Util::output_spectrum_mo_librpa(velocity_mo,
                                           outdir + "velocity_matrix",
                                           nk,
                                           nspin_tmp,
                                           nbands,
                                           nbasis,
                                           nbands,
                                           *this->p_kv);
    }
    ModuleBase::timer::end("RPA_LRI", "out_velocity");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_librpa_basis_v1(const UnitCell& ucell,
                                            Exx_LRI<double>* exx_lri,
                                            const std::string& aux_filename,
                                            const std::string& legacy_filename)
{
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_librpa_basis_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> type_naux(static_cast<std::size_t>(ucell.ntype), 0);
    std::vector<int> type_nw(static_cast<std::size_t>(ucell.ntype), 0);
    int total_wfc = 0;
    int total_aux = 0;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        type_naux[static_cast<std::size_t>(it)] = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(it);
        if (type_naux[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type auxiliary size.");
        }
        type_nw[static_cast<std::size_t>(it)] = ucell.atoms[it].nw;
        if (type_nw[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type wave-function size.");
        }
        total_wfc += type_nw[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
        total_aux += type_naux[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
    }

    const auto wfc_l_nchi = RpaLriDetail::collect_wfc_l_nchi(ucell);
    const auto aux_l_nchi = RpaLriDetail::collect_abfs_l_nchi(exx_lri->abfs);
    const std::string wfc_filename = outdir + "basis_wfc_out";
    const std::string aux_output_filename = outdir + aux_filename;
    const std::string legacy_output_filename = outdir + legacy_filename;
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_nw, wfc_l_nchi, wfc_filename);
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_naux, aux_l_nchi, aux_output_filename);

    std::ofstream ofs(legacy_output_filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + legacy_output_filename);
    }
    ofs << std::setw(10) << ucell.ntype
        << std::setw(10) << total_wfc
        << std::setw(10) << total_aux
        << "    fallback" << std::endl;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        ofs << std::setw(10) << it + 1
            << std::setw(10) << type_nw[static_cast<std::size_t>(it)]
            << std::setw(10) << type_naux[static_cast<std::size_t>(it)]
            << std::endl;
    }
    const auto write_legacy_shells = [&ofs](const std::vector<std::vector<int>>& shell_counts_by_type)
    {
        for (std::size_t itype = 0; itype < shell_counts_by_type.size(); ++itype)
        {
            const auto& shell_counts = shell_counts_by_type[itype];
            int nshell = 0;
            for (const int count : shell_counts)
            {
                nshell += count;
            }
            ofs << std::setw(10) << itype + 1
                << std::setw(10) << nshell
                << std::endl;
            for (std::size_t l = 0; l < shell_counts.size(); ++l)
            {
                for (int iradial = 0; iradial < shell_counts[l]; ++iradial)
                {
                    ofs << std::setw(10) << l << std::endl;
                }
            }
        }
    };
    write_legacy_shells(wfc_l_nchi);
    write_legacy_shells(aux_l_nchi);
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k_v1(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                         std::string filename,
                                         Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k_v1");
    ModuleBase::timer::tick("RPA_LRI", "out_coulomb_k_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    const auto atom_naux = this->collect_atom_naux_(ucell, exx_lri);
    const int all_mu = RpaLriDetail::sum_int_vector(atom_naux);
    const int nks_tot = PARAM.inp.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const std::size_t natoms = static_cast<std::size_t>(ucell.nat);

    struct V1Block
    {
        int pair_index = 0;
        int I = 0;
        int J = 0;
        std::int64_t offset = 0;
        RI::Tensor<std::complex<double>> tensor;
    };

    for (int ik = 0; ik != nks_tot; ++ik)
    {
        std::vector<V1Block> blocks;
        for (auto& Ip: Vs)
        {
            const int I = static_cast<int>(Ip.first);
            const int mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                const int J = static_cast<int>(JPp.first.first);
                if (J < I)
                {
                    continue;
                }
                if (!RpaLriDetail::has_valid_matrix_shape(JPp.second))
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const auto R = JPp.first.second;
                const double arg = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                    * ModuleBase::TWO_PI;
                const std::complex<double> kphase = std::complex<double>(std::cos(arg), std::sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                const int J = static_cast<int>(vq_Jp.first);
                auto& vq_J = vq_Jp.second;
                const int nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[J]);
                if (static_cast<int>(vq_J.shape[0]) != mu_num || static_cast<int>(vq_J.shape[1]) != nu_num)
                {
                    throw std::runtime_error("LibRPA v1 Coulomb output encountered an inconsistent tensor shape.");
                }
                V1Block block;
                block.pair_index = static_cast<int>(RpaLriDetail::coulomb_atom_pair_index(
                    static_cast<std::size_t>(I), static_cast<std::size_t>(J), natoms));
                block.I = I;
                block.J = J;
                block.tensor = std::move(vq_J);
                blocks.push_back(std::move(block));
            }
        }

        std::sort(blocks.begin(), blocks.end(), [](const V1Block& lhs, const V1Block& rhs) {
            return lhs.pair_index < rhs.pair_index;
        });
        if (blocks.empty())
        {
            continue;
        }
        for (std::size_t ib = 1; ib < blocks.size(); ++ib)
        {
            if (blocks[ib - 1].pair_index == blocks[ib].pair_index)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output found duplicate atom-pair blocks on one MPI rank.");
            }
        }

        const std::int64_t nblocks = static_cast<std::int64_t>(blocks.size());
        const std::int64_t header_size = static_cast<std::int64_t>(6 * sizeof(std::int32_t))
            + static_cast<std::int64_t>(ucell.nat * sizeof(std::int32_t))
            + nblocks * static_cast<std::int64_t>(sizeof(std::int32_t) + sizeof(std::int64_t));
        std::int64_t offset = header_size;
        for (auto& block: blocks)
        {
            block.offset = offset;
            const unsigned long long values = RpaLriDetail::checked_mul_u64(
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.I)]),
                static_cast<unsigned long long>(atom_naux[static_cast<std::size_t>(block.J)]),
                "LibRPA v1 Coulomb block size");
            const unsigned long long bytes = RpaLriDetail::checked_mul_u64(
                values,
                static_cast<unsigned long long>(sizeof(std::complex<double>)),
                "LibRPA v1 Coulomb block size");
            offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Coulomb block size");
        }

        std::stringstream ss;
        ss << filename << ik + 1 << "_rank" << GlobalV::MY_RANK << ".dat";
        const std::string out_name = outdir + ss.str();
        std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + out_name);
        }

        const std::int32_t marker = RpaLriDetail::LIBRPA_COULOMB_V1_MARKER;
        const std::int32_t iq = ik + 1;
        const std::int32_t naux = all_mu;
        const std::int32_t value_flag = RpaLriDetail::LIBRPA_COULOMB_V1_COMPLEX_FLAG;
        const std::int32_t natom = ucell.nat;
        const std::int32_t nblock_i32 = static_cast<std::int32_t>(nblocks);
        RpaLriDetail::write_scalar(ofs, marker, out_name);
        RpaLriDetail::write_scalar(ofs, iq, out_name);
        RpaLriDetail::write_scalar(ofs, naux, out_name);
        RpaLriDetail::write_scalar(ofs, value_flag, out_name);
        RpaLriDetail::write_scalar(ofs, natom, out_name);
        RpaLriDetail::write_scalar(ofs, nblock_i32, out_name);
        for (const int atom_aux: atom_naux)
        {
            const std::int32_t atom_aux_i32 = atom_aux;
            RpaLriDetail::write_scalar(ofs, atom_aux_i32, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::int32_t pair_index = block.pair_index;
            RpaLriDetail::write_scalar(ofs, pair_index, out_name);
            RpaLriDetail::write_scalar(ofs, block.offset, out_name);
        }
        for (const auto& block: blocks)
        {
            const std::size_t nvalues = block.tensor.get_shape_all();
            if (nvalues == 0)
            {
                throw std::runtime_error("LibRPA v1 Coulomb output encountered an empty tensor payload.");
            }
            RpaLriDetail::checked_write(
                ofs,
                block.tensor.ptr(),
                nvalues * sizeof(std::complex<double>),
                out_name);
        }
        ofs.close();
    }

    ModuleBase::timer::tick("RPA_LRI", "out_coulomb_k_v1");
}


// template<typename Tdata>
// void RPA_LRI<T, Tdata>::init(const MPI_Comm &mpi_comm_in)
// {
// 	if(this->info == this->exx.info)
// 	{
// 		this->lcaos = this->exx.lcaos;
// 		this->abfs = this->exx.abfs;
// 		this->abfs_ccp = this->exx.abfs_ccp;

// 		exx_lri_rpa.cv = std::move(this->exx.cv);
// 	}
// 	else
// 	{
// 		this->lcaos = ...
// 		this->abfs = ...
// 		this->abfs_ccp = ...

// 		exx_lri_rpa.cv.set_orbitals(
// 			this->lcaos, this->abfs, this->abfs_ccp,
// 			this->info.kmesh_times, this->info.ccp_rmesh_times );
// 	}

// }



// template<typename Tdata>
// void RPA_LRI<T, Tdata>::cal_rpa_ions()
// {
// 	// this->rpa_lri.set_parallel(this->mpi_comm, atoms_pos, latvec, period);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Vws = std::move(this->exx.cv.Vws);

// 	const std::array<Tcell,Ndim> period_Vs =
// LRI_CV_Tools::cal_latvec_range<Tcell>(1+this->info.ccp_rmesh_times); const
// std::pair<std::vector<TA>,
// std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>> 		list_As_Vs
// = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs,
// 2, false);

// 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>
// 		Vs = exx_lri_rpa.cv.cal_Vs(
// 			list_As_Vs.first, list_As_Vs.second[0],
// 			{{"writable_Vws",true}});

// 	// Vs[iat0][{iat1,cell1}]	distributed across processes by (iat0,iat1); each process holds all cell1
// 	Vqs = FFT(Vs);
// 	out_Vs(Vqs);

// 	if(this->info == this->exx.info)
// 		exx_lri_rpa.cv.Cws = std::move(this->exx.cv.Cws);

// 	const std::array<Tcell,Ndim> period_Cs =
// LRI_CV_Tools::cal_latvec_range<Tcell>(2); 	const std::pair<std::vector<TA>,
// std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>> 		list_As_Cs
// = RI::Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms,
// period_Cs, 2, false);

// 	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,
// std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>> 		Cs_dCs =
// exx_lri_rpa.cv.cal_Cs_dCs( 			list_As_Cs.first, list_As_Cs.second[0],
// 			{{"cal_dC",false},
// 			 {"writable_Cws",true}, {"writable_dCws",true},
// {"writable_Vws",false},
// {"writable_dVws",false}}); 	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs
// = std::get<0>(Cs_dCs);

// 	out_Cs(Cs);

// 	// rpa_lri.set_Cs(Cs);
// }

#endif
