#ifndef LAT_R_CSR_H
#define LAT_R_CSR_H

#include "hs_sparse_io.h"

#include <fstream>

class Parallel_Blkcyclic; // Defined in source/source_base/parallel_blkcyc.h

namespace ModuleIO
{
    template <typename T>
    void save_lat_r(std::ofstream& ofs,
        const SparseRBlock<T>& XR,
        const Parallel_Blkcyclic& layout,
        const SparseWriteOptions& options);
}

#endif
