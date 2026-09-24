//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

// Born-von Karmen supercell utilities.
// This header is free of LibRI dependencies so that modules built without
// LibRI (e.g. module_lr) can use it.

#ifndef RI_UTIL_BVK_H
#define RI_UTIL_BVK_H

#include "source_cell/klist.h"

#include <array>
#include <vector>

namespace RI_Util
{
    inline std::array<int,3>
    get_Born_vonKarmen_period(const K_Vectors &kv)
    {
        return std::array<int,3>{kv.nmp[0], kv.nmp[1], kv.nmp[2]};
    }

    // Fold a cell index into [-period/2, period/2) to match LibRI's
    // Array_Operator::operator%: (c % period + 3 * period / 2) % period - period / 2
    template<typename Tcell>
    inline Tcell fold_cell_centered(const Tcell c, const Tcell period)
    {
        return (c % period + 3 * period / 2) % period - period / 2;
    }

    template<typename Tcell>
    std::vector<std::array<Tcell,1>>
    get_Born_von_Karmen_cells( const std::array<Tcell,1> &Born_von_Karman_period )
    {
        std::vector<std::array<Tcell,1>> Born_von_Karman_cells;
        for( Tcell c=0; c<Born_von_Karman_period[0]; ++c )
            Born_von_Karman_cells.emplace_back(
                std::array<Tcell,1>{ fold_cell_centered(c, Born_von_Karman_period[0]) });
        return Born_von_Karman_cells;
    }

    template<typename Tcell, size_t Ndim>
    std::vector<std::array<Tcell,Ndim>>
    get_Born_von_Karmen_cells( const std::array<Tcell,Ndim> &Born_von_Karman_period )
    {
        std::array<Tcell,Ndim-1> sub_Born_von_Karman_period;
        for(size_t i=0; i<Ndim-1; ++i)
            sub_Born_von_Karman_period[i] = Born_von_Karman_period[i];

        std::vector<std::array<Tcell,Ndim>> Born_von_Karman_cells;
        for( const std::array<Tcell,Ndim-1> &sub_cell : get_Born_von_Karmen_cells(sub_Born_von_Karman_period) )
            for( Tcell c=0; c<Born_von_Karman_period.back(); ++c )
            {
                std::array<Tcell,Ndim> cell;
                for(size_t i=0; i<Ndim-1; ++i)
                    cell[i] = sub_cell[i];
                cell.back() = fold_cell_centered(c, Born_von_Karman_period.back());
                Born_von_Karman_cells.emplace_back(std::move(cell));
            }
        return Born_von_Karman_cells;
    }

    /* example for Ndim=3:
    template<typename Tcell, size_t Ndim>
    std::vector<std::array<Tcell,Ndim>>
    get_Born_von_Karmen_cells( const std::array<Tcell,Ndim> &Born_von_Karman_period )
    {
        std::vector<std::array<Tcell,Ndim>> Born_von_Karman_cells;
        for( int ix=0; ix<Born_von_Karman_period[0]; ++ix )
            for( int iy=0; iy<Born_von_Karman_period[1]; ++iy )
                for( int iz=0; iz<Born_von_Karman_period[2]; ++iz )
                    Born_von_Karman_cells.push_back( std::array<Tcell,Ndim>{ix,iy,iz} );
        return Born_von_Karman_cells;
    }
    */
}

#endif
