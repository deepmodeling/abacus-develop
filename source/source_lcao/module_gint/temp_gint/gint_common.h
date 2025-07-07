#pragma once
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_gint/temp_gint/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    void compose_hr_gint(HContainer<double>& hr_gint);
    // for nspin=4 case
    void compose_hr_gint(const std::vector<HContainer<double>>& hr_gint_part,
                         HContainer<std::complex<double>>& hr_gint_full);

    template <typename T>
    void transfer_hr_gint_to_hR(const HContainer<T>& hr_gint, HContainer<T>& hR);

    template<typename T>
    void transfer_dm_2d_to_gint(
        const GintInfo& gint_info,
        std::vector<HContainer<T>*> dm,
        std::vector<HContainer<T>>& dm_gint);

    template<typename T>
    void wfc_2d_to_gint(const T* wfc_2d, int nbands, int nlocal, const Parallel_Orbitals& pv, T* wfc_grid, const GintInfo& gint_info);
}
