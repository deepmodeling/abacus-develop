#include "biggrid_info.h"
#include "gint_helper.h"
#include "gint_type.h"

namespace ModuleGint
{

BigGridInfo::BigGridInfo(
    Vec3d biggrid_vec1,
    Vec3d biggrid_vec2,
    Vec3d biggrid_vec3,
    int nmx, int nmy, int nmz)
    : biggrid_vec1_(biggrid_vec1),
      biggrid_vec2_(biggrid_vec2),
      biggrid_vec3_(biggrid_vec3),
      nmx_(nmx), nmy_(nmy), nmz_(nmz), nmxyz_(nmx*nmy*nmz)
    {
        // initialize the biggrid_latvec0_
        biggrid_latvec0_.e11 = biggrid_vec1_.x;
        biggrid_latvec0_.e12 = biggrid_vec1_.y;
        biggrid_latvec0_.e13 = biggrid_vec1_.z;

        biggrid_latvec0_.e21 = biggrid_vec2_.x;
        biggrid_latvec0_.e22 = biggrid_vec2_.y;
        biggrid_latvec0_.e23 = biggrid_vec2_.z;

        biggrid_latvec0_.e31 = biggrid_vec3_.x;
        biggrid_latvec0_.e32 = biggrid_vec3_.y;
        biggrid_latvec0_.e33 = biggrid_vec3_.z;

        // initialize the GT matrix
        biggrid_GT_ = biggrid_latvec0_.Inverse();

        // initialize the meshgrid_info_
        meshgrid_info_ = std::make_shared<MeshGridInfo>(
            biggrid_vec1_ / static_cast<double>(nmx),
            biggrid_vec2_ / static_cast<double>(nmy),
            biggrid_vec3_ / static_cast<double>(nmz));
        
        // initialize the meshgrid_coords_
        meshgrid_coords_.resize(nmxyz_);
        for(int index_1d = 0; index_1d < nmxyz_; index_1d++)
        {
            meshgrid_coords_[index_1d] = 
                meshgrid_info_->get_cartesian_coord(mgrid_idx_1Dto3D(index_1d));
        }
    }

    Vec3i BigGridInfo::max_ext_bgrid_num(double r) const
    {
        const double g1 = sqrt(biggrid_GT_.e11 * biggrid_GT_.e11
            + biggrid_GT_.e21 * biggrid_GT_.e21
            + biggrid_GT_.e31 * biggrid_GT_.e31);
        const double g2 = sqrt(biggrid_GT_.e12 * biggrid_GT_.e12
            + biggrid_GT_.e22 * biggrid_GT_.e22
            + biggrid_GT_.e32 * biggrid_GT_.e32);
        const double g3 = sqrt(biggrid_GT_.e13 * biggrid_GT_.e13
            + biggrid_GT_.e23 * biggrid_GT_.e23
            + biggrid_GT_.e33 * biggrid_GT_.e33);
        int ext_x = static_cast<int>(r * g1) + 1;
        int ext_y = static_cast<int>(r * g2) + 1;
        int ext_z = static_cast<int>(r * g3) + 1;
        return Vec3i(ext_x, ext_y, ext_z);
    }

} // namespace ModuleGint