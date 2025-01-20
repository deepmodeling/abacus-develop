#pragma once
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_gint/temp_gint/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hRGint(std::shared_ptr<HContainer<double>> hRGint);
    // for nspin=4 case
    void compose_hRGint(std::vector<std::shared_ptr<HContainer<double>>> hRGint_part,
        std::shared_ptr<HContainer<std::complex<double>>> hRGint_full);

    template <typename T>
    void transfer_hRGint_to_hR(std::shared_ptr<const HContainer<T>> hRGint, HContainer<T>* hR);

    void transfer_DM_to_DMGint(
        std::shared_ptr<const GintInfo> gint_info,
        std::vector<HContainer<double>*> DM,
        std::vector<std::shared_ptr<HContainer<double>>> DMRGint);

}
