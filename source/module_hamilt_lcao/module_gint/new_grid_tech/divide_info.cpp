#include "divide_info.h"

namespace ModuleGint
{

DivideInfo::DivideInfo(
    int startidx_bx_old, int startidx_by_old, int startidx_bz_old,
    int nbx_old, int nby_old, int nbz_old,
    std::shared_ptr<const UnitCellInfo> unitcell_info, bool is_redevided)
    : startidx_bx_old_(startidx_bx_old), startidx_by_old_(startidx_by_old), startidx_bz_old_(startidx_bz_old),
      nbx_old_(nbx_old), nby_old_(nby_old), nbz_old_(nbz_old),
      startidx_bx_new_(startidx_bx_old), startidx_by_new_(startidx_by_old), startidx_bz_new_(startidx_bz_old),
      nbx_new_(nbx_old), nby_new_(nby_old), nbz_new_(nbz_old),
      unitcell_info_(unitcell_info), is_redivided_(is_redevided)
    {
        if(!is_redivided_)
        {
            localcell_info_ = std::make_shared<LocalCellInfo>(startidx_bx_new_, startidx_by_new_, startidx_bz_new_,
                                                            nbx_new_, nby_new_, nbz_new_, unitcell_info_);
        }
        // TODO: "implement the redivide function";
    }

}