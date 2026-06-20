#ifndef SINGLE_R_IO_H
#define SINGLE_R_IO_H

#include "write_HS_sparse.h"

#include <fstream>

namespace ModuleIO
{
    template <typename T>
    void output_single_R(std::ofstream& ofs,
        const SparseRBlock<T>& XR,
        const Parallel_Orbitals& pv,
        const SparseWriteOptions& options);
}

#endif
