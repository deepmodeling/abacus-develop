#include "psi_init_nao_random.h"

template <typename T>
void psi_init_nao_random<T>::init_psig(T* psig, const int& ik)
{
    double rm = this->mixing_coef_;
    psi_init_nao<T>::init_psig(psig, ik);
    const int npol = this->npol_;
    const int nbasis = this->pw_wfc_->npwk_max * npol;
    psi::Psi<T> psi_random(1, this->nbands_start_, nbasis, nbasis, true);
    psi_random.fix_k(0);
    this->random_t(psi_random.get_pointer(), 0, this->nbands_start_, ik, 0);
    for (int iband = 0; iband < this->nbands_start_; iband++)
    {
        for (int ibasis = 0; ibasis < nbasis; ibasis++)
        {
            psig[iband * nbasis + ibasis] *= (T(1.0) + Real(rm) * psi_random(iband, ibasis));
        }
    }
}

template class psi_init_nao_random<std::complex<double>>;
template class psi_init_nao_random<std::complex<float>>;
// gamma point calculation
template class psi_init_nao_random<double>;
template class psi_init_nao_random<float>;