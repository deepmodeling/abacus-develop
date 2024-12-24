#pragma once
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hRGint(hamilt::HContainer<double>* hR);

    template <typename T>
    void transfer_hRGint_to_hR(const hamilt::HContainer<T>* hRGint, hamilt::HContainer<T>* hR);

}
