#ifndef RECORD_ADJ_H
#define RECORD_ADJ_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

//---------------------------------------------------
// FUNCTION: record the adjacent atoms for each atom
//---------------------------------------------------
class Record_adj
{
  private:
    bool info_modified = false;

  public:
    Record_adj();
    ~Record_adj();

    /**
     * Build the complete per-center adjacency records used by LCAO consumers.
     *
     * With an MPI Parallel_Orbitals communicator, spatial owners perform the
     * neighbor filtering and exchange compact records before each rank rebuilds
     * the complete global-iat view. The 2D orbital distribution is then used to
     * update pv.nlocstart, pv.nlocdim and pv.nnr.
     *
     * mpi_grid may point to a prebuilt local-center Grid with the same build
     * radius, physical query radius and boundary mode as grid_d. Passing it
     * avoids rebuilding MPI domain data inside Record_adj. If omitted in an
     * MPI run, a compatible temporary distributed Grid is built collectively.
     */
    void for_2d(const UnitCell& ucell,
                const Grid_Driver& grid_d,
                Parallel_Orbitals& pv,
                bool gamma_only,
                const std::vector<double>& orb_cutoff,
                const Grid_Driver* mpi_grid = nullptr);

    void delete_grid();

    int na_proc = 0;
    int* na_each = nullptr;

    //--------------------------------------------
    // record sparse atom index in for_grid();
    // Map iat(dense atom index) to sparse atom index
    // Mainly removing the index dependency for OpenMP parallel loop
    //
    // Meaning:
    // 1. if iat2ca[iat] > 0, it contains the sparse atom index
    // 2. if iat2ca[iat] < 0, the sparse atom index of iat does not exist
    //
    // Usage:
    // 1. iat2ca[iat] > 0 ? na_each[iat2ca[iat]] : 0
    // 2. iat2ca[iat] > 0 ? info[iat2ca[iat]] : nullptr
    //--------------------------------------------
    int* iat2ca = nullptr;

    //------------------------------------------------
    // info will identify each atom in each unitcell.
    //------------------------------------------------
    int*** info = nullptr;

  private:
};

#endif
