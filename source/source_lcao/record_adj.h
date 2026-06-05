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

    //--------------------------------------------
    // This will record the orbitals according to
    // HPSEPS's 2D block division.
    //--------------------------------------------
    void for_2d(const UnitCell& ucell,
                const Grid_Driver& grid_d,
                Parallel_Orbitals& pv,
                bool gamma_only,
                const std::vector<double>& orb_cutoff);


    void delete_grid();

    int na_proc=0;
    int* na_each=nullptr;

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
    int* iat2ca=nullptr;

    //------------------------------------------------
    // info will identify each atom in each unitcell.
    //------------------------------------------------
    int*** info=nullptr;

  private:
};

#endif
