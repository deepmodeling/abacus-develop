#include "biggrid_info.h"
#include "gint_helper.h"

namespace Gint
{

BigGridInfo::BigGridInfo(
    const Vector3& biggrid_vec1,
    const Vector3& biggrid_vec2,
    const Vector3& biggrid_vec3,
    const int nmx, const int nmy, const int nmz)
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
                this->get_meshgrid_coord(index_1d);
        }
    }

} // namespace Gint