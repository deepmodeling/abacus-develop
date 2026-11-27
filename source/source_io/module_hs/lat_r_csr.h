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
     * @param layout Block-cyclic distribution used for the MPI reduction.
     * @param options Sparse output options (threshold, precision, reduction).
     */
    template <typename T>
    void save_lat_r_csr(std::ofstream& ofs,
        const SparseRBlock<T>& XR,
        const Parallel_Blkcyclic& layout,
        const SparseWriteOptions& options);

    /**
     * @brief Write an all-zero R block in the H(R)/S(R)/r(R) CSR text layout.
     * @details Callers that already know a block is globally empty (for
     *          example from an allreduced nonzero count) use this instead of
     *          save_lat_r_csr to skip its row-wise MPI reductions. The emitted
     *          bytes are identical to save_lat_r_csr on an empty block: an
     *          empty SparseMatrix printed on the rank that owns the stream.
     * @param ofs Output text stream.
     * @param nlocal Number of localized basis functions.
     * @param options Sparse output options (precision and reduce flag).
     */
    template <typename T>
    void save_lat_r_empty_csr(std::ofstream& ofs,
        int nlocal,
        const SparseWriteOptions& options);
}

#endif
