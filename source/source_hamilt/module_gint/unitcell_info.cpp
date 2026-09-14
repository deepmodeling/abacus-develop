#include "unitcell_info.h"
#include "gint_helper.h"

namespace ModuleGint
{
    
UnitCellInfo::UnitCellInfo(
    const Vec3d& unitcell_vec1,
    const Vec3d& unitcell_vec2,
    const Vec3d& unitcell_vec3,
    int nbx, int nby, int nbz,
    int nmx, int nmy, int nmz)
    : unitcell_vec1_(unitcell_vec1),
        unitcell_vec2_(unitcell_vec2),
        unitcell_vec3_(unitcell_vec3),
        nbx_(nbx), nby_(nby), nbz_(nbz), nbxyz_(nbx*nby*nbz),
        nmx_(nmx), nmy_(nmy), nmz_(nmz), nmxyz_(nmx*nmy*nmz)
    {
        // initialize the biggrid_info_
        biggrid_info_ = std::make_shared<BigGridInfo>(
            unitcell_vec1_ / static_cast<double>(nbx),
            unitcell_vec2_ / static_cast<double>(nby),
            unitcell_vec3_ / static_cast<double>(nbz),
            nmx/nbx, nmy/nby, nmz/nbz);
        
        meshgrid_info_ = biggrid_info_->get_mgrid_info();
    }
    
}