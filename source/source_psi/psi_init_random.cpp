#include "psi_init_random.h"

template <typename T>
void psi_init_random<T>::initialize(const Structure_Factor* sf,
                                               const ModulePW::PW_Basis_K* pw_wfc,
                                               const UnitCell* p_ucell,
                                               const K_Vectors* p_kv_in,
                                               const int& random_seed,
                                               const pseudopot_cell_vnl* p_pspot_nl,
                                               const int& rank)
{
    psi_initializer<T>::initialize(sf, pw_wfc, p_ucell, p_kv_in, random_seed, p_pspot_nl, rank);
}

template <typename T>
void psi_init_random<T>::init_psig(T* psig, const int& ik)
{
    this->random_t(psig, 0, this->nbands_start_, ik, 0);
}

template class psi_init_random<std::complex<double>>;
template class psi_init_random<std::complex<float>>;
// gamma point calculation
template class psi_init_random<double>;
template class psi_init_random<float>;