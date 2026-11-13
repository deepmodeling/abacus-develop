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

    /**
     * @brief Write the CSR-format description banner used by H(R)/S(R)/r(R)
     *        sparse text output.
     * @param ofs Output text stream.
     */
    void write_csr_banner(std::ofstream& ofs);

    /**
     * @brief Write one R block in the H(R)/S(R)/r(R) CSR text layout.
     * @details The block is reduced across ranks exactly like save_lat_r, then
     *          printed through SparseMatrix so that values wrap every six
     *          entries and an empty block still carries its three CSR
     *          sections. Only the root rank opens the output stream, so
     *          `ofs.is_open()` selects the rank that owns the assembled
     *          matrix and performs the printing; every rank still takes part
     *          in the collective reduction.
     * @param ofs Output text stream.
     * @param XR Sparse block of the selected matrix component.
     * @param pv Parallel orbital distribution used for the MPI reduction.
     * @param options Sparse output options (threshold, precision, reduction).
     */
    template <typename T>
    void save_lat_r_csr(std::ofstream& ofs,
        const SparseRBlock<T>& XR,
        const Parallel_Orbitals& pv,
        const SparseWriteOptions& options);
}

#endif
