#include "source_base/module_parallel/para_band_output.h"

#include "source_base/global_function.h"
#include "source_base/module_parallel/para_mpi_func.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace Parallel
{

ParaBandOutput::ParaBandOutput(const int local_nbands, const int global_nbands, const ParaWorld& band_world)
    : band_world_(band_world), global_nbands_(global_nbands), local_nbands_(local_nbands)
{
    if (local_nbands < 0 || global_nbands < 0)
    {
        ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput", "band counts cannot be negative");
    }
    if (!this->band_world_.valid())
    {
        ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput", "band domain must be valid");
    }
    const int band_groups = this->band_world_.size();
    this->band_group_ = this->band_world_.rank();
    this->band_counts_.resize(band_groups);
    // The domain connects corresponding k-point and PW ranks across band groups.
    int band_count = local_nbands;
    allgather_int(band_count, this->band_counts_.data(), this->band_world_);

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
            ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput", "local band counts do not match global nbands");
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
}

bool ParaBandOutput::bands_are_replicated() const
{
    return this->bands_are_replicated_;
}

int ParaBandOutput::global_nbands() const
{
    return this->global_nbands_;
}

int ParaBandOutput::local_nbands() const
{
    return this->local_nbands_;
}

int ParaBandOutput::band_group() const
{
    return this->band_group_;
}

int ParaBandOutput::local_offset() const
{
    return this->band_offsets_[this->band_group_];
}

int ParaBandOutput::owner_group(const int global_band) const
{
    if (global_band < 0 || global_band >= this->global_nbands_)
    {
        ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput", "global band index is out of range");
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
    ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput", "global band has no owner");
    return 0;
}

int ParaBandOutput::local_index(const int global_band) const
{
    const int owner = this->owner_group(global_band);
    return global_band - this->band_offsets_[owner];
}

void ParaBandOutput::bcast_band(const int global_band, double* data, const int count) const
{
    const int owner = this->owner_group(global_band);
    bcast_double(data, count, this->band_world_, owner);
}

void ParaBandOutput::bcast_band(const int global_band, std::complex<double>* data, const int count) const
{
    const int owner = this->owner_group(global_band);
    bcast_complex(data, count, this->band_world_, owner);
}

ModuleBase::matrix ParaBandOutput::gather_matrix(const ModuleBase::matrix& local_matrix) const
{
    // Every rank must agree to proceed before any rank returns or gathers rows.
    int invalid_columns = local_matrix.nc != this->local_nbands_ ? 1 : 0;
    reduce_max(invalid_columns, this->band_world_);
    if (invalid_columns != 0)
    {
        ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput::gather_matrix", "matrix columns do not match the band layout");
    }
    if (this->bands_are_replicated_)
    {
        // Keep local replicas unchanged, including the empty-band case.
        return local_matrix;
    }
    const int band_groups = this->band_world_.size();
    std::vector<int> row_counts(band_groups);
    int local_rows = local_matrix.nr;
    allgather_int(local_rows, row_counts.data(), this->band_world_);

    const bool rows_match = std::all_of(row_counts.begin(), row_counts.end(), [&local_matrix](const int rows) { return rows == local_matrix.nr; });
    if (!rows_match)
    {
        ModuleBase::WARNING_QUIT("Parallel::ParaBandOutput::gather_matrix", "band groups have inconsistent row counts");
    }

    ModuleBase::matrix global_matrix(local_matrix.nr, this->global_nbands_, false);
    // Rows represent local k-points and columns represent contiguous band shards.
    for (int ik = 0; ik < local_matrix.nr; ++ik)
    {
        const double* local_row = local_matrix.nc > 0 ? local_matrix.c + ik * local_matrix.nc : nullptr;
        allgatherv_double(local_row,
                          local_matrix.nc,
                          global_matrix.c + ik * this->global_nbands_,
                          this->band_counts_.data(),
                          this->band_offsets_.data(),
                          this->band_world_);
    }
    return global_matrix;
}

} // namespace Parallel
