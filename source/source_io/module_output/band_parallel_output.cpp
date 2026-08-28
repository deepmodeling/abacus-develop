#include "band_parallel_output.h"

#include "source_base/global_function.h"
#include "source_base/parallel_comm.h"

#include <algorithm>
#include <numeric>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

ModuleIO::BandParallelLayout::BandParallelLayout(const int local_nbands, const int global_nbands)
    : global_nbands_(global_nbands), local_nbands_(local_nbands)
{
    if (local_nbands < 0 || global_nbands < 0)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::BandParallelLayout", "band counts cannot be negative");
    }
#ifndef __MPI
    // A serial build has no communicator that could supply missing band shards.
    if (local_nbands != global_nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::BandParallelLayout", "distributed bands require an MPI build");
    }
    this->band_counts_.push_back(local_nbands);
    this->band_offsets_.push_back(0);
#else
    int band_groups = 0;
    MPI_Comm_size(BP_WORLD, &band_groups);
    MPI_Comm_rank(BP_WORLD, &this->band_group_);
    this->band_counts_.resize(band_groups);
    // BP_WORLD connects ranks with the same k-point and plane-wave slab but different
    // band groups, so its rank order is also the band-group ownership order.
    MPI_Allgather(&local_nbands, 1, MPI_INT, this->band_counts_.data(), 1, MPI_INT, BP_WORLD);

    // SDFT may replicate all deterministic bands in every group. BPCG instead stores
    // complementary contiguous shards whose total must equal the global band count.
    this->bands_are_replicated_ = std::all_of(this->band_counts_.begin(), this->band_counts_.end(), [global_nbands](const int bands) {
        return bands == global_nbands;
    });
    if (!this->bands_are_replicated_)
    {
        const int gathered_nbands = std::accumulate(this->band_counts_.begin(), this->band_counts_.end(), 0);
        if (gathered_nbands != global_nbands)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::BandParallelLayout", "local band counts do not match global nbands");
        }
    }

    this->band_offsets_.resize(band_groups, 0);
    if (!this->bands_are_replicated_)
    {
        // BPCG assigns shards in band-group order, making the prefix sum both the
        // global starting offset and the basis for global-to-local mapping.
        for (int group = 1; group < band_groups; ++group)
        {
            this->band_offsets_[group] = this->band_offsets_[group - 1] + this->band_counts_[group - 1];
        }
    }
#endif
}

bool ModuleIO::BandParallelLayout::bands_are_replicated() const
{
    return this->bands_are_replicated_;
}

int ModuleIO::BandParallelLayout::global_nbands() const
{
    return this->global_nbands_;
}

int ModuleIO::BandParallelLayout::local_nbands() const
{
    return this->local_nbands_;
}

int ModuleIO::BandParallelLayout::band_group() const
{
    return this->band_group_;
}

int ModuleIO::BandParallelLayout::local_offset() const
{
    return this->band_offsets_[this->band_group_];
}

int ModuleIO::BandParallelLayout::owner_group(const int global_band) const
{
    if (global_band < 0 || global_band >= this->global_nbands_)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::BandParallelLayout", "global band index is out of range");
    }
    if (this->bands_are_replicated_)
    {
        // Designate group 0 as the output owner to avoid redundant work and writes.
        return 0;
    }
    // Offsets delimit half-open contiguous ownership ranges [offset, offset + count).
    for (int group = 0; group < static_cast<int>(this->band_counts_.size()); ++group)
    {
        if (global_band < this->band_offsets_[group] + this->band_counts_[group])
        {
            return group;
        }
    }
    ModuleBase::WARNING_QUIT("ModuleIO::BandParallelLayout", "global band has no owner");
    return 0;
}

int ModuleIO::BandParallelLayout::local_index(const int global_band) const
{
    const int owner = this->owner_group(global_band);
    return global_band - this->band_offsets_[owner];
}

ModuleBase::matrix ModuleIO::gather_band_matrix(const ModuleBase::matrix& local_matrix, const int global_nbands)
{
    const BandParallelLayout layout(local_matrix.nc, global_nbands);
    if (layout.bands_are_replicated_)
    {
        // The local matrix is already complete; no BP_WORLD data movement is needed.
        return local_matrix;
    }
#ifdef __MPI
    const int band_groups = static_cast<int>(layout.band_counts_.size());
    std::vector<int> row_counts(band_groups);
    MPI_Allgather(&local_matrix.nr, 1, MPI_INT, row_counts.data(), 1, MPI_INT, BP_WORLD);

    const bool rows_match
        = std::all_of(row_counts.begin(), row_counts.end(), [&local_matrix](const int rows) { return rows == local_matrix.nr; });
    if (!rows_match)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::gather_band_matrix", "band groups have inconsistent row counts");
    }

    ModuleBase::matrix global_matrix(local_matrix.nr, global_nbands, false);
    // Rows represent local k-points and columns represent band shards. Gather each row
    // into band-group order so every BP_WORLD rank receives the same global matrix.
    for (int ik = 0; ik < local_matrix.nr; ++ik)
    {
        const double* local_row = local_matrix.nc > 0 ? local_matrix.c + ik * local_matrix.nc : nullptr;
        MPI_Allgatherv(local_row,
                       local_matrix.nc,
                       MPI_DOUBLE,
                       global_matrix.c + ik * global_nbands,
                       layout.band_counts_.data(),
                       layout.band_offsets_.data(),
                       MPI_DOUBLE,
                       BP_WORLD);
    }
    return global_matrix;
#else
    return local_matrix;
#endif
}
