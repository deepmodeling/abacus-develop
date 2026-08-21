#include "band_parallel_output.h"

#include "source_base/global_function.h"
#include "source_base/parallel_comm.h"

#include <numeric>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

ModuleBase::matrix ModuleIO::gather_band_matrix(const ModuleBase::matrix& local_matrix, const int global_nbands)
{
    assert(global_nbands > 0);
#ifndef __MPI
    if (local_matrix.nc != global_nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::gather_band_matrix", "distributed bands require an MPI build");
    }
    return local_matrix;
#else
    int band_groups = 0;
    MPI_Comm_size(BP_WORLD, &band_groups);
    if (band_groups == 1)
    {
        if (local_matrix.nc != global_nbands)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::gather_band_matrix", "local band count does not match global nbands");
        }
        return local_matrix;
    }

    std::vector<int> band_counts(band_groups);
    MPI_Allgather(&local_matrix.nc, 1, MPI_INT, band_counts.data(), 1, MPI_INT, BP_WORLD);
    const int gathered_nbands = std::accumulate(band_counts.begin(), band_counts.end(), 0);
    if (gathered_nbands != global_nbands)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::gather_band_matrix", "local band counts do not match global nbands");
    }

    std::vector<int> band_offsets(band_groups, 0);
    for (int group = 1; group < band_groups; ++group)
    {
        band_offsets[group] = band_offsets[group - 1] + band_counts[group - 1];
    }

    ModuleBase::matrix global_matrix(local_matrix.nr, global_nbands, false);
    for (int ik = 0; ik < local_matrix.nr; ++ik)
    {
        const double* local_row = local_matrix.nc > 0 ? local_matrix.c + ik * local_matrix.nc : nullptr;
        MPI_Allgatherv(local_row,
                       local_matrix.nc,
                       MPI_DOUBLE,
                       global_matrix.c + ik * global_nbands,
                       band_counts.data(),
                       band_offsets.data(),
                       MPI_DOUBLE,
                       BP_WORLD);
    }
    return global_matrix;
#endif
}
