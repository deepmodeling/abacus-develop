#pragma once
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_gint/new_grid_tech/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hRGint(hamilt::HContainer<double>* hR);

    template <typename T>
    void transfer_hRGint_to_hR(const hamilt::HContainer<T>* hRGint, hamilt::HContainer<T>* hR);

    void transfer_DM_to_DMGint(
        std::shared_ptr<const GintInfo> gint_info,
        std::vector<hamilt::HContainer<double>*> DM,
        std::vector<std::shared_ptr<hamilt::HContainer<double>>> DMRGint);

}
