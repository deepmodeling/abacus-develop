#pragma once

#include <memory>
#include <vector>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/matrix.h"
#include "gint.h"
#include "gint_info.h"

namespace ModuleGint
{
class Gint_fvl_meta : public Gint
{
    public:
    Gint_fvl_meta(
        const int nspin,
        const std::vector<const double*>& vr_eff,
        const std::vector<const double*>& vofk,
        const std::vector<HContainer<double>*>& dm_vec,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl,
        ModuleBase::matrix* svl)
        : nspin_(nspin), vr_eff_(vr_eff), vofk_(vofk), dm_vec_(dm_vec),
          isforce_(isforce), isstress_(isstress), fvl_(fvl), svl_(svl),
          dr3_(gint_info_->get_mgrid_volume()) {};

    void cal_gint() override;

    private:
    void init_dm_gint_();

    void cal_fvl_svl_();

    // input
    const int nspin_;
    std::vector<const double*> vr_eff_;
    std::vector<const double*> vofk_;
    std::vector<HContainer<double>*> dm_vec_;
    const bool isforce_;
    const bool isstress_;

    // output
    ModuleBase::matrix* fvl_;
    ModuleBase::matrix* svl_;

    // intermediate variables
    std::vector<std::shared_ptr<HContainer<double>>> dm_gint_vec_;

    double dr3_;
};

} // namespace ModuleGint