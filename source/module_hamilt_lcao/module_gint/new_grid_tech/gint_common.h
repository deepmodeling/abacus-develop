#pragma once
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_gint/new_grid_tech/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hRGint(std::shared_ptr<hamilt::HContainer<double>> hRGint);
    // for nspin=4 case
    void compose_hRGint(std::vector<std::shared_ptr<hamilt::HContainer<double>>> hRGint_part,
        std::shared_ptr<hamilt::HContainer<std::complex<double>>> hRGint_full);

    template <typename T>
    void transfer_hRGint_to_hR(std::shared_ptr<const hamilt::HContainer<T>> hRGint, hamilt::HContainer<T>* hR);

    void transfer_DM_to_DMGint(
        std::shared_ptr<const GintInfo> gint_info,
        std::vector<hamilt::HContainer<double>*> DM,
        std::vector<std::shared_ptr<hamilt::HContainer<double>>> DMRGint);

}
