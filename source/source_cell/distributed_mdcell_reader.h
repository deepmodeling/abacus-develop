#ifndef DISTRIBUTED_MDCELL_READER_H
#define DISTRIBUTED_MDCELL_READER_H

#include "source_cell/md_cell.h"

#ifdef __MPI
#include <mpi.h>
#include <string>

class DistributedMdCellReader
{
public:
    static MdCell read_lj_stru(const std::string& stru_file,
                               MPI_Comm comm,
                               double cutoff_bohr,
                               double skin_bohr = 0.0);
};

#endif

#endif
