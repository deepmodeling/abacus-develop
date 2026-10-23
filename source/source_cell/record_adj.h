#ifndef RECORD_ADJ_H
#define RECORD_ADJ_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/unitcell.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

#include <array>
#include <vector>

//---------------------------------------------------
// FUNCTION: record the adjacent atoms for each atom
//---------------------------------------------------
class Record_adj
{
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
                const int npol,
                const std::vector<double>& orb_cutoff);


    void delete_grid();

  private:
    // (1) count adjacent atoms for each atom and,
    // for multi-k, accumulate nlocdim / nlocstart / nnr of pv.
    void count_adjacent(const UnitCell& ucell,
                        const Grid_Driver& grid_d,
                        Parallel_Orbitals& pv,
                        bool gamma_only,
                        const int npol,
                        const std::vector<double>& orb_cutoff);

    // allocate info[na_proc][na_each[i]][5]
    void allocate_info();

    // fill info with (Rx, Ry, Rz, T, I) of each adjacent atom.
    void fill_info(const UnitCell& ucell,
                   const Grid_Driver& grid_d,
                   const std::vector<double>& orb_cutoff);

  public:
    int na_proc=0;
    std::vector<int> na_each;

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
    std::vector<int> iat2ca;

    //------------------------------------------------
    // info identifies each adjacent atom in each
    // unitcell. All adjacent records are stored flat:
    // the records of atom iat occupy
    // info[info_offset[iat], info_offset[iat]+na_each[iat]).
    // Each record holds (Rx, Ry, Rz, T, I).
    //------------------------------------------------
    std::vector<std::array<int, 5>> info;
    std::vector<int> info_offset;

    // Access the (Rx, Ry, Rz, T, I) record of the cb-th
    // adjacent atom of atom iat.
    const std::array<int, 5>& get_info(const int iat, const int cb) const
    {
        return info[info_offset[iat] + cb];
    }
};

#endif
