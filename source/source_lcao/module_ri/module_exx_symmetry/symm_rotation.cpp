#include "symm_rotation.h"

#include "symm_rotation_r_detail.h"

#include <algorithm>

namespace
{
// Convert one cached ComplexMatrix block to the RI::Tensor the real-space
// rotation code needs. File-local: only get_rotmat_Slm_tensor calls it.
RI::Tensor<std::complex<double>> rotmat_Slm_to_tensor(const ModuleBase::ComplexMatrix& block)
    {
        RI::Tensor<std::complex<double>> t({ static_cast<size_t>(block.nr), static_cast<size_t>(block.nc) });
        for (int i = 0; i < block.nr; ++i) {
            for (int j = 0; j < block.nc; ++j) {
                t(i, j) = block(i, j);
        } }
        return t;
    }
} // namespace

namespace ModuleSymmetry
{
    void Symmetry_rotation::set_Cs_rotation(const std::vector<std::vector<int>>& abfs_l_nchi)
    {
        this->reduce_Cs_ = true;
        this->abfs_l_nchi_ = abfs_l_nchi;
        for (const std::vector<int>& abfs_T : abfs_l_nchi) { this->abfs_Lmax_ = std::max(this->abfs_Lmax_, static_cast<int>(abfs_T.size()) - 1); }
    }

    /// Cached (isym, l) lookup into rotmat_Slm_tensor_: rebuilt only when rotmat_Slm_version_
    /// (bumped by cal_rotmat_Slm, i.e. once per ion step) changes, instead of reconverting the
    /// same block on every set_rotation_matrix/set_rotation_matrix_abf call (once per atom pair).
    const RI::Tensor<std::complex<double>>& Symmetry_rotation::get_rotmat_Slm_tensor(const int isym, const int l)const
    {
        if (this->rotmat_Slm_tensor_version_ != this->rotmat_Slm_version_)
        {
            this->rotmat_Slm_tensor_.resize(this->rotmat_Slm_.size());
            for (size_t is = 0; is < this->rotmat_Slm_.size(); ++is)
            {
                this->rotmat_Slm_tensor_[is].resize(this->rotmat_Slm_[is].size());
                for (size_t il = 0; il < this->rotmat_Slm_[is].size(); ++il) {
                    this->rotmat_Slm_tensor_[is][il] = rotmat_Slm_to_tensor(this->rotmat_Slm_[is][il]);
            } }
            this->rotmat_Slm_tensor_version_ = this->rotmat_Slm_version_;
        }
        return this->rotmat_Slm_tensor_[isym][l];
    }

}
