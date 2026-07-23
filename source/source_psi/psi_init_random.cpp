#include "psi_init_random.h"
#include <vector>
#include "source_io/module_parameter/parameter.h"

template <typename T>
void psi_init_random<T>::initialize(const Structure_Factor* sf,
                                               const ModulePW::PW_Basis_K* pw_wfc,
                                               const UnitCell* p_ucell,
                                               const std::vector<int>& ik2iktot,
                                               const int& nkstot,
                                               const int& random_seed,
                                               const pseudopot_cell_vnl* p_pspot_nl,
                                               const int& rank,
                                               const int& npol,
                                               const int& nbands)
{
    psi_initializer<T>::initialize(sf, pw_wfc, p_ucell, ik2iktot, nkstot, random_seed, p_pspot_nl, rank, npol, nbands);
    this->ixy2is_.clear();
    this->ixy2is_.resize(this->pw_wfc_->fftnxy);
    this->pw_wfc_->getfftixy2is(this->ixy2is_.data());
    this->nbands_start_ = nbands;
    this->nbands_complem_ = 0;
}

template <typename T>
void psi_init_random<T>::init_psig(T* psig, const int& ik)
{
    this->random_t(psig, 0, this->nbands_start_, ik);
}

template class psi_init_random<std::complex<double>>;
template class psi_init_random<std::complex<float>>;
// gamma point calculation
template class psi_init_random<double>;
template class psi_init_random<float>;
