#pragma once
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hr_gint(std::shared_ptr<HContainer<double>> hr_gint);
    // for nspin=4 case
    void compose_hr_gint(std::vector<std::shared_ptr<HContainer<double>>> hr_gint_part,
        std::shared_ptr<HContainer<std::complex<double>>> hr_gint_full);

    template <typename T>
    void transfer_hr_gint_to_hR(std::shared_ptr<const HContainer<T>> hr_gint, HContainer<T>* hR);

    void transfer_dm_2d_to_gint(
        std::shared_ptr<const GintInfo> gint_info,
        std::vector<HContainer<double>*> dm,
        std::vector<std::shared_ptr<HContainer<double>>> dm_gint);

}
