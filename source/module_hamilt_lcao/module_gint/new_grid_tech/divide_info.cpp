#include "devide_info.h"

namespace Gint
{

DivideInfo::DivideInfo(
    int startidx_bx_old, int startidx_by_old, int startidx_bz_old,
    int nbx_old, int nby_old, int nbz_old,
    std::shared_ptr<const UnitCellInfo> unitcell_info, bool is_redevided)
    : startidx_bx_old(startidx_bx_old), startidx_by_old(startidx_by_old), startidx_bz_old(startidx_bz_old),
      nbx_old(nbx_old), nby_old(nby_old), nbz_old(nbz_old),
      startidx_bx_new(startidx_bx_old), startidx_by_new(startidx_by_old), startidx_bz_new(startidx_bz_old),
      nbx_new(nbx_old), nby_new(nby_old), nbz_new(nbz_old),
      unitcell_info(unitcell_info), is_redevided(is_redevided)
    {
        if(!is_redvided)
        {
            localcell_info = std::make_shared<LocalCellInfo>(startidx_bx_new_, startidx_by_new_, startidx_bz_new_,
                                                            nbx_new_, nby_new_, nbz_new_, unitcell_info_);
        }
        // TODO: "implement the redivide function";
    }

}