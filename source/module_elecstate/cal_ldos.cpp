#include "elecstate_lcao.h"
#include "elecstate_pw.h"

namespace elecstate
{
template <typename T, typename Device>
void ElecStatePW<T, Device>::cal_ldos(const psi::Psi<T, Device>& psi, std::vector<double>& ldos) const
{
    // energy range for ldos (efermi as reference)
    const double emin = PARAM.inp.stm_bias < 0 ? PARAM.inp.stm_bias : 0;
    const double emax = PARAM.inp.stm_bias > 0 ? PARAM.inp.stm_bias : 0;

    std::vector<T> wfcr(this->basis->nrxx);

    for (int ik = 0; ik < this->klist->get_nks(); ++ik)
    {
        psi.fix_k(ik);
        double efermi = this->eferm.get_efval(this->klist->isk[ik]);
        int nbands = psi.get_nbands();

        for (int ib = 0; ib < nbands; ib++)
        {
            this->basis->recip2real(&psi(ib, 0), wfcr.data(), ik);
            double eigenval = (ekb(ik, ib) - efermi) * ModuleBase::Ry_to_eV;
            if (eigenval >= emin && eigenval <= emax)
            {
                for (int ir = 0; ir < basis->nrxx; ir++)
                    ldos[ir] += this->klist->wk[ik] * norm(wfcr[ir]);
            }
        }
    }
}

template class ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

// lcao multi-k case
template <>
void ElecStateLCAO<std::complex<double>>::cal_ldos(const psi::Psi<std::complex<double>>& psi,
                                                   std::vector<double>& ldos) const
{
}

// lcao Gamma_only case
template <>
void ElecStateLCAO<double>::cal_ldos(const psi::Psi<double>& psi, std::vector<double>& ldos) const
{
}

template class ElecStateLCAO<double>;               // Gamma_only case
template class ElecStateLCAO<std::complex<double>>; // multi-k case
} // namespace elecstate
