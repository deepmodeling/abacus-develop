#pragma once
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_gint/gint_info.h"

namespace ModuleGint
{
    // fill the lower triangle matrix with the upper triangle matrix
    template<typename T>
    void compose_hr_gint(HContainer<T>& hr_gint);

    template<typename Tout, typename Tin>
    void cast_hcontainer_values(const HContainer<Tin>& src, HContainer<Tout>& dst);

    template<typename Tout, typename Tin>
    HContainer<Tout> make_cast_hcontainer(const HContainer<Tin>& src);
    

    template <typename T>
    void transfer_hr_gint_to_hR(const HContainer<T>& hr_gint, HContainer<T>& hR);
    // for nspin=4 case
    void merge_hr_part_to_hR(const std::vector<hamilt::HContainer<double>>& hr_gint_tmp ,
                         hamilt::HContainer<std::complex<double>>* hR,
                         const GintInfo& gint_info);

    template<typename TGint, typename TDM>
    void transfer_dm_2d_to_gint(
        const GintInfo& gint_info,
        const std::vector<HContainer<TDM>*>& dm,
        std::vector<HContainer<TGint>>& dm_gint);

    template<typename T>
    void wfc_2d_to_gint(const T* wfc_2d, int nbands, int nlocal, const Parallel_Orbitals& pv, T* wfc_grid, const GintInfo& gint_info);
}
