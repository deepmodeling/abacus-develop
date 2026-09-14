#pragma once

#include "biggrid_info.h"
#include "unitcell_info.h"
#include "localcell_info.h"

namespace ModuleGint
{

class DivideInfo
{
    public:
    // constructor
    DivideInfo(
        int startidx_bx_old, int startidx_by_old, int startidx_bz_old,
        int nbx_old, int nby_old, int nbz_old,
        std::shared_ptr<const UnitCellInfo> unitcell_info, bool is_redivided = false);
    
    // getter functions
    std::shared_ptr<const LocalCellInfo> get_localcell_info() const { return localcell_info_; }
    bool get_is_redivided() const { return is_redivided_; }
    
    private:
    // if the grid is redivided, is_redeiided_ is true
    bool is_redivided_;
    
    // the old start index of the local cell
    int startidx_bx_old_;
    int startidx_by_old_;
    int startidx_bz_old_;

    // the old number of big grids in the local cell
    int nbx_old_;
    int nby_old_;
    int nbz_old_;

    // the new start index of the local cell
    int startidx_bx_new_;
    int startidx_by_new_;
    int startidx_bz_new_;

    // the new number of big grids in the local cell
    int nbx_new_;
    int nby_new_;
    int nbz_new_;

    // the unitcell info
    std::shared_ptr<const UnitCellInfo> unitcell_info_;

    // the localcell info
    std::shared_ptr<const LocalCellInfo> localcell_info_;
};

}