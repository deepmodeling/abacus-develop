#include "sternheimer_siab_mpi.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace module_ri
{
namespace sternheimer_siab
{
namespace
{

using Complex = std::complex<double>;

std::vector<ReferenceRow> sorted_rows(std::vector<ReferenceRow> rows)
{
    const auto key = [](const ReferenceRow& row) {
        return std::make_tuple(row.occupied_state, row.auxiliary_channel, row.frequency_index);
    };
    std::sort(rows.begin(), rows.end(), [&key](const ReferenceRow& left, const ReferenceRow& right) {
        return key(left) < key(right);
    });
    for (std::size_t index = 1; index < rows.size(); ++index)
    {
        if (key(rows[index - 1]) == key(rows[index]))
        {
            throw std::invalid_argument("Sternheimer SIAB gathered reference rows contain a duplicate key.");
        }
    }
    return rows;
}

bool valid_local_rows(const std::vector<ReferenceRow>& rows, const std::size_t nprimitive)
{
    for (const ReferenceRow& row: rows)
    {
        if (row.q.size() != nprimitive || row.occupied_state < 0 || row.auxiliary_channel < 0 || row.frequency_index < 0
            || !std::isfinite(row.frequency_ha) || !std::isfinite(row.occupation)
            || !std::isfinite(row.frequency_weight) || !std::isfinite(row.norm))
        {
            return false;
        }
        for (const Complex& value: row.q)
        {
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                return false;
            }
        }
    }
    return true;
}

std::size_t row_width(const std::size_t nprimitive)
{
    if (nprimitive > (std::numeric_limits<std::size_t>::max() - 7) / 2)
    {
        throw std::overflow_error("Sternheimer SIAB reference row width overflows size_t.");
    }
    return 7 + 2 * nprimitive;
}

std::vector<double> pack_rows(const std::vector<ReferenceRow>& rows, const std::size_t nprimitive)
{
    const std::size_t width = row_width(nprimitive);
    std::vector<double> packed;
    packed.reserve(rows.size() * width);
    for (const ReferenceRow& row: rows)
    {
        packed.push_back(static_cast<double>(row.occupied_state));
        packed.push_back(static_cast<double>(row.auxiliary_channel));
        packed.push_back(static_cast<double>(row.frequency_index));
        packed.push_back(row.frequency_ha);
        packed.push_back(row.occupation);
        packed.push_back(row.frequency_weight);
        packed.push_back(row.norm);
        for (const Complex& value: row.q)
        {
            packed.push_back(value.real());
            packed.push_back(value.imag());
        }
    }
    return packed;
}

int unpack_index(const double value)
{
    if (!std::isfinite(value) || value < 0.0 || value > static_cast<double>(std::numeric_limits<int>::max())
        || value != std::floor(value))
    {
        throw std::runtime_error("Sternheimer SIAB gathered reference row contains an invalid integer index.");
    }
    return static_cast<int>(value);
}

std::vector<ReferenceRow> unpack_rows(const std::vector<double>& packed, const std::size_t nprimitive)
{
    const std::size_t width = row_width(nprimitive);
    if (packed.size() % width != 0)
    {
        throw std::runtime_error("Sternheimer SIAB gathered reference row payload has an invalid size.");
    }
    std::vector<ReferenceRow> rows(packed.size() / width);
    for (std::size_t irow = 0; irow != rows.size(); ++irow)
    {
        const std::size_t offset = irow * width;
        ReferenceRow& row = rows[irow];
        row.occupied_state = unpack_index(packed[offset]);
        row.auxiliary_channel = unpack_index(packed[offset + 1]);
        row.frequency_index = unpack_index(packed[offset + 2]);
        row.frequency_ha = packed[offset + 3];
        row.occupation = packed[offset + 4];
        row.frequency_weight = packed[offset + 5];
        row.norm = packed[offset + 6];
        row.q.resize(nprimitive);
        for (std::size_t ie = 0; ie != nprimitive; ++ie)
        {
            row.q[ie] = Complex(packed[offset + 7 + 2 * ie], packed[offset + 8 + 2 * ie]);
        }
    }
    return sorted_rows(std::move(rows));
}

} // namespace

std::vector<std::vector<Complex>> assemble_full_primitive_grids(const std::vector<PrimitiveSlab>& slabs,
                                                                const int nxy,
                                                                const int nz)
{
    if (nxy <= 0 || nz <= 0 || slabs.empty())
    {
        throw std::invalid_argument("Sternheimer SIAB primitive assembly requires positive grid dimensions and slabs.");
    }
    const std::size_t nprimitive = slabs.front().values.size();
    if (nprimitive == 0)
    {
        throw std::invalid_argument("Sternheimer SIAB primitive assembly requires at least one primitive.");
    }

    std::vector<int> plane_coverage(static_cast<std::size_t>(nz), 0);
    std::vector<std::vector<Complex>> full(
        nprimitive,
        std::vector<Complex>(static_cast<std::size_t>(nxy) * static_cast<std::size_t>(nz), Complex(0.0, 0.0)));
    for (const PrimitiveSlab& slab: slabs)
    {
        if (slab.startz < 0 || slab.nplane < 0 || slab.startz + slab.nplane > nz || slab.values.size() != nprimitive)
        {
            throw std::invalid_argument("Sternheimer SIAB primitive slab metadata is inconsistent.");
        }
        const std::size_t local_size = static_cast<std::size_t>(nxy) * static_cast<std::size_t>(slab.nplane);
        for (const auto& primitive: slab.values)
        {
            if (primitive.size() != local_size)
            {
                throw std::invalid_argument("Sternheimer SIAB primitive slab payload size is inconsistent.");
            }
        }
        for (int iz = 0; iz != slab.nplane; ++iz)
        {
            ++plane_coverage[static_cast<std::size_t>(slab.startz + iz)];
        }
        for (std::size_t ie = 0; ie != nprimitive; ++ie)
        {
            for (int ixy = 0; ixy != nxy; ++ixy)
            {
                for (int iz = 0; iz != slab.nplane; ++iz)
                {
                    full[ie][static_cast<std::size_t>(ixy) * static_cast<std::size_t>(nz)
                             + static_cast<std::size_t>(slab.startz + iz)]
                        = slab.values[ie][static_cast<std::size_t>(ixy) * static_cast<std::size_t>(slab.nplane)
                                          + static_cast<std::size_t>(iz)];
                }
            }
        }
    }
    if (!std::all_of(plane_coverage.begin(), plane_coverage.end(), [](const int count) { return count == 1; }))
    {
        throw std::invalid_argument("Sternheimer SIAB primitive slabs must cover every PW z plane exactly once.");
    }
    return full;
}

#ifdef __MPI
std::vector<std::vector<Complex>> allgather_full_primitive_grids(
    const std::vector<std::vector<Complex>>& local_primitives,
    const int nxy,
    const int nz,
    const int startz,
    const int nplane,
    MPI_Comm communicator)
{
    int rank_count = 0;
    MPI_Comm_size(communicator, &rank_count);
    const bool local_valid
        = nxy > 0 && nz > 0 && startz >= 0 && nplane >= 0 && startz + nplane <= nz && !local_primitives.empty()
          && std::all_of(local_primitives.begin(), local_primitives.end(), [nxy, nplane](const auto& values) {
                 return values.size() == static_cast<std::size_t>(nxy) * static_cast<std::size_t>(nplane);
             });
    int valid = local_valid ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, communicator);
    if (valid == 0)
    {
        throw std::invalid_argument("Sternheimer SIAB local primitive slab is invalid on at least one MPI rank.");
    }

    const int local_nprimitive = static_cast<int>(local_primitives.size());
    int min_nprimitive = local_nprimitive;
    int max_nprimitive = local_nprimitive;
    MPI_Allreduce(MPI_IN_PLACE, &min_nprimitive, 1, MPI_INT, MPI_MIN, communicator);
    MPI_Allreduce(MPI_IN_PLACE, &max_nprimitive, 1, MPI_INT, MPI_MAX, communicator);
    if (min_nprimitive != max_nprimitive)
    {
        throw std::invalid_argument("Sternheimer SIAB primitive counts differ between MPI ranks.");
    }

    const std::size_t local_grid_size = static_cast<std::size_t>(nxy) * static_cast<std::size_t>(nplane);
    const std::size_t local_payload_size = local_grid_size * local_primitives.size();
    if (local_payload_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Sternheimer SIAB local primitive payload exceeds MPI int count capacity.");
    }
    std::vector<Complex> local_payload;
    local_payload.reserve(local_payload_size);
    for (const auto& primitive: local_primitives)
    {
        local_payload.insert(local_payload.end(), primitive.begin(), primitive.end());
    }

    std::vector<int> startz_by_rank(static_cast<std::size_t>(rank_count));
    std::vector<int> nplane_by_rank(static_cast<std::size_t>(rank_count));
    std::vector<int> counts(static_cast<std::size_t>(rank_count));
    const int local_count = static_cast<int>(local_payload_size);
    MPI_Allgather(&startz, 1, MPI_INT, startz_by_rank.data(), 1, MPI_INT, communicator);
    MPI_Allgather(&nplane, 1, MPI_INT, nplane_by_rank.data(), 1, MPI_INT, communicator);
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator);

    std::vector<int> displacements(static_cast<std::size_t>(rank_count), 0);
    std::size_t total_count = 0;
    for (int rank = 0; rank != rank_count; ++rank)
    {
        if (total_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::overflow_error("Sternheimer SIAB primitive allgather exceeds MPI int displacement capacity.");
        }
        displacements[static_cast<std::size_t>(rank)] = static_cast<int>(total_count);
        total_count += static_cast<std::size_t>(counts[static_cast<std::size_t>(rank)]);
    }
    std::vector<Complex> gathered(total_count);
    MPI_Allgatherv(local_payload.data(),
                   local_count,
                   MPI_DOUBLE_COMPLEX,
                   gathered.data(),
                   counts.data(),
                   displacements.data(),
                   MPI_DOUBLE_COMPLEX,
                   communicator);

    std::vector<PrimitiveSlab> slabs(static_cast<std::size_t>(rank_count));
    for (int rank = 0; rank != rank_count; ++rank)
    {
        PrimitiveSlab& slab = slabs[static_cast<std::size_t>(rank)];
        slab.startz = startz_by_rank[static_cast<std::size_t>(rank)];
        slab.nplane = nplane_by_rank[static_cast<std::size_t>(rank)];
        const std::size_t rank_grid_size = static_cast<std::size_t>(nxy) * static_cast<std::size_t>(slab.nplane);
        slab.values.resize(static_cast<std::size_t>(local_nprimitive));
        for (int ie = 0; ie != local_nprimitive; ++ie)
        {
            const std::size_t begin = static_cast<std::size_t>(displacements[static_cast<std::size_t>(rank)])
                                      + static_cast<std::size_t>(ie) * rank_grid_size;
            slab.values[static_cast<std::size_t>(ie)].assign(gathered.begin() + begin,
                                                             gathered.begin() + begin + rank_grid_size);
        }
    }
    return assemble_full_primitive_grids(slabs, nxy, nz);
}

std::vector<ReferenceRow> gather_reference_rows_to_root(const std::vector<ReferenceRow>& local_rows,
                                                        const std::size_t nprimitive,
                                                        const int root,
                                                        MPI_Comm communicator)
{
    int rank = 0;
    int rank_count = 0;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &rank_count);
    if (root < 0 || root >= rank_count)
    {
        throw std::invalid_argument("Sternheimer SIAB row gather root is out of range.");
    }
    int valid = valid_local_rows(local_rows, nprimitive) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, communicator);
    if (valid == 0)
    {
        throw std::invalid_argument("Sternheimer SIAB local reference rows are invalid on at least one MPI rank.");
    }

    const std::vector<double> local_packed = pack_rows(local_rows, nprimitive);
    if (local_packed.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Sternheimer SIAB local row payload exceeds MPI int count capacity.");
    }
    const int local_count = static_cast<int>(local_packed.size());
    std::vector<int> counts(rank == root ? static_cast<std::size_t>(rank_count) : 0);
    MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, root, communicator);

    std::vector<int> displacements;
    std::vector<double> gathered;
    if (rank == root)
    {
        displacements.assign(static_cast<std::size_t>(rank_count), 0);
        std::size_t total_count = 0;
        for (int source_rank = 0; source_rank != rank_count; ++source_rank)
        {
            if (total_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            {
                throw std::overflow_error("Sternheimer SIAB row gather exceeds MPI int displacement capacity.");
            }
            displacements[static_cast<std::size_t>(source_rank)] = static_cast<int>(total_count);
            total_count += static_cast<std::size_t>(counts[static_cast<std::size_t>(source_rank)]);
        }
        gathered.resize(total_count);
    }
    MPI_Gatherv(local_packed.data(),
                local_count,
                MPI_DOUBLE,
                gathered.data(),
                counts.data(),
                displacements.data(),
                MPI_DOUBLE,
                root,
                communicator);
    return rank == root ? unpack_rows(gathered, nprimitive) : std::vector<ReferenceRow>();
}
#else
std::vector<std::vector<Complex>> allgather_full_primitive_grids(
    const std::vector<std::vector<Complex>>& local_primitives,
    const int nxy,
    const int nz,
    const int startz,
    const int nplane)
{
    PrimitiveSlab slab;
    slab.startz = startz;
    slab.nplane = nplane;
    slab.values = local_primitives;
    return assemble_full_primitive_grids({slab}, nxy, nz);
}

std::vector<ReferenceRow> gather_reference_rows_to_root(const std::vector<ReferenceRow>& local_rows,
                                                        const std::size_t nprimitive,
                                                        const int root)
{
    if (root != 0 || !valid_local_rows(local_rows, nprimitive))
    {
        throw std::invalid_argument("Sternheimer SIAB serial reference row gather input is invalid.");
    }
    return sorted_rows(local_rows);
}
#endif

} // namespace sternheimer_siab
} // namespace module_ri
