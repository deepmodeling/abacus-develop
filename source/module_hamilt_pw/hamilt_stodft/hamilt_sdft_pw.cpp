#include "hamilt_sdft_pw.h"
#include "module_base/timer.h"

namespace hamilt
{

template <typename T, typename Device>
HamiltSdftPW<T, Device>::HamiltSdftPW(elecstate::Potential* pot_in,
                                      ModulePW::PW_Basis_K* wfc_basis,
                                      K_Vectors* p_kv,
                                      const int& npol,
                                      double* emin_in,
                                      double* emax_in)
    : HamiltPW<T, Device>(pot_in, wfc_basis, p_kv), ngk(p_kv->ngk)
{
    this->classname = "HamiltSdftPW";
    this->npwk_max = wfc_basis->npwk_max;
    this->npol = npol;
    this->emin = emin_in;
    this->emax = emax_in;
}

template <typename T, typename Device>
void HamiltSdftPW<T, Device>::hPsi(const T* psi_in, T* hpsi, const int& nbands)
{
    auto call_act = [&, this](const Operator<T, Device>* op) -> void {
        op->act(nbands, this->npwk_max, this->npol, psi_in, hpsi, this->ngk[op->get_ik()]);
    };

    ModuleBase::timer::tick("HamiltSdftPW", "hPsi");
    ModuleBase::GlobalFunc::ZEROS(hpsi, nbands * this->npwk_max * this->npol);
    call_act(this->ops);
    Operator<T, Device>* node((Operator<T, Device>*)this->ops->next_op);
    while (node != nullptr)
    {
        call_act(node);
        node = (Operator<T, Device>*)(node->next_op);
    }
    ModuleBase::timer::tick("HamiltSdftPW", "hPsi");

    return;
}

template <typename T, typename Device>
void HamiltSdftPW<T, Device>::hPsi_norm(const T* psi_in, T* hpsi_norm, const int& nbands)
{
    ModuleBase::timer::tick("HamiltSdftPW", "hPsi_norm");

    this->hPsi(psi_in, hpsi_norm, nbands);

    const int ik = this->ops->get_ik();
    const int npwk_max = this->npwk_max;
    const int npwk = this->ngk[ik];
    using Real = typename GetTypeReal<T>::type;
    const Real emin = *this->emin;
    const Real emax = *this->emax;
    const Real Ebar = (emin + emax) / 2;
    const Real DeltaE = (emax - emin) / 2;
    for (int ib = 0; ib < nbands; ++ib)
    {
        for (int ig = 0; ig < npwk; ++ig)
        {
            hpsi_norm[ib * npwk_max + ig]
                = (hpsi_norm[ib * npwk_max + ig] - Ebar * psi_in[ib * npwk_max + ig]) / DeltaE;
        }
    }
    ModuleBase::timer::tick("HamiltSdftPW", "hPsi_norm");
}

template class HamiltSdftPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HamiltSdftPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HamiltSdftPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HamiltSdftPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt