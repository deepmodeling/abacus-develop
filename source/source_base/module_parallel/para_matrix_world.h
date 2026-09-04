#ifndef PARA_MATRIX_WORLD_H
#define PARA_MATRIX_WORLD_H

#include <cstddef>
#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief matrix parallel domain: 2D block-cyclic distribution.
 *
 * Self-contained wrapper for the matrix-level parallel topology.
 * Replaces the process-grid part of Parallel_2D (dim0/dim1/coord)
 * with a cleaner interface. The actual ScaLAPACK descriptor and
 * BLACS context management stay in Parallel_2D; this class only
 * holds the process grid dimensions and coordinates.
 *
 * Tests only need this header.
 */
class ParaMatrixWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial matrix domain (1x1 process grid).
     */
    ParaMatrixWorld();

#ifdef __MPI
    /**
     * @brief Construct a matrix domain on an existing communicator.
     *
     * The process grid is computed automatically: dim0 = largest divisor
     * of nproc with dim0 >= dim1 (square-ish), dim1 = nproc / dim0.
     *
     * @param[in] comm  matrix communicator (e.g. DIAG_WORLD or MPI_COMM_WORLD)
     */
    ParaMatrixWorld(const MPI_Comm& comm);
#endif

    /// Process grid row count.
    int dim0() const { return dim0_; }

    /// Process grid column count.
    int dim1() const { return dim1_; }

    /// Row coordinate of this process in the grid.
    int coord_row() const { return coord_row_; }

    /// Column coordinate of this process in the grid.
    int coord_col() const { return coord_col_; }

private:
    int dim0_ = 1;
    int dim1_ = 1;
    int coord_row_ = 0;
    int coord_col_ = 0;

    void compute_proc_grid();
};

} // namespace Parallel

#endif // PARA_MATRIX_WORLD_H
