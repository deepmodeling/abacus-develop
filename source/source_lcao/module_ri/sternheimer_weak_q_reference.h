#ifndef STERNHEIMER_WEAK_Q_REFERENCE_H
#define STERNHEIMER_WEAK_Q_REFERENCE_H

#include "source_lcao/module_ri/sternheimer_abacus_st_smoke.h"

#include <cstdint>
#include <string>
#include <vector>

namespace ModuleRI
{

struct SternheimerWeakQReferenceDimensions
{
    std::uint64_t fine_grid_size = 0;
    std::uint64_t full_record_count = 0;
    std::uint64_t orbital_count = 0;
};

struct SternheimerWeakQReferenceLimits
{
    std::uint64_t max_file_bytes = UINT64_C(1073741824);
    // Cumulative decoded element storage, including nested vector objects;
    // not a bound on allocator bookkeeping or the caller's existing copies.
    std::uint64_t max_allocation_bytes = UINT64_C(1073741824);
    std::uint64_t max_bands_per_record = UINT64_C(1048576);
};

struct SternheimerWeakQReference
{
    std::vector<double> fine_potential;
    std::vector<SternheimerLCAOOccupiedKPoint> full_records;
};

// contract_hash is an opaque, nonempty byte string (at most 65536 bytes).
// The caller must hash its exact physical contract: geometry, grids, PP/NAO/
// orbital ordering, spin, units, etc. No normalization or SCF-value hashing is
// performed here. Dimensions must be positive and match exactly on both ends.
//
// Schema 1: little-endian fixed-width integers and IEEE binary64, all actual
// occupied-k-point fields, followed by 64 lowercase ASCII SHA256 hex bytes.
// The checksum covers the entire preceding file, including contract and schema.
// This detects corruption, not malicious substitution by a party able to reseal.
// Bump the schema and update both codecs when the record definition changes.
//
// POSIX publication uses a same-directory temporary file, fsync, and link: an
// existing destination (including a dangling symlink) is never overwritten.
// All validation/I/O failures throw std::runtime_error; no partial reference
// is returned or published. Allocator exhaustion may still throw std::bad_alloc.
void write_sternheimer_weak_q_reference(
    const std::string& path,
    const std::string& contract_hash,
    const SternheimerWeakQReferenceDimensions& dimensions,
    const std::vector<double>& fine_potential,
    const std::vector<SternheimerLCAOOccupiedKPoint>& full_records,
    const SternheimerWeakQReferenceLimits& limits = {});

// After successful validation, the caller must replace BOTH its fine-potential
// copy and its full records before computing the old reference hash, e.g.:
//   auto reference = read_sternheimer_weak_q_reference(path, contract, dims);
//   fine_potential.swap(reference.fine_potential);
//   full_records.swap(reference.full_records);
// The two swaps are nonthrowing with these standard-allocator vector types.
SternheimerWeakQReference read_sternheimer_weak_q_reference(
    const std::string& path,
    const std::string& contract_hash,
    const SternheimerWeakQReferenceDimensions& dimensions,
    const SternheimerWeakQReferenceLimits& limits = {});

} // namespace ModuleRI

#endif
