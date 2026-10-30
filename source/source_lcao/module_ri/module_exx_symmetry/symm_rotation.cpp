#include "symm_rotation.h"

#include <algorithm>

namespace ModuleSymmetry
{
    void Symmetry_rotation::set_Cs_rotation(const std::vector<std::vector<int>>& abfs_l_nchi)
    {
        this->reduce_Cs_ = true;
        this->abfs_l_nchi_ = abfs_l_nchi;
        for (auto& abfs_T : abfs_l_nchi) { this->abfs_Lmax_ = std::max(this->abfs_Lmax_, static_cast<int>(abfs_T.size()) - 1); }
    }

}
