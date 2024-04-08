#ifndef SINGLE_R_IO_H
#define SINGLE_R_IO_H

#include "module_basis/module_ao/parallel_orbitals.h"

namespace ModuleIO
{
    template <typename T>
    void output_single_R(std::ofstream& ofs,
        const std::map<size_t, std::map<size_t, T>>& XR,
        const double& sparse_threshold,
        const bool& binary,
        const Parallel_Orbitals& pv,
        const bool& reduce = true);
}

#endif
