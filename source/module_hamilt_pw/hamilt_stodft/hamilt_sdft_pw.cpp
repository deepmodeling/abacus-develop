#include "hamilt_sdft_pw.h"
#include "module_base/timer.h"
#include "kernels/hpsi_norm_op.h"

namespace hamilt
{

template <typename T, typename Device>
HamiltSdftPW<T, Device>::HamiltSdftPW(elecstate::Potential* pot_in,
                                      ModulePW::PW_Basis_K* wfc_basis,
                                      K_Vectors* p_kv,
                                      pseudopot_cell_vnl* nlpp,
                                      const UnitCell* ucell,
                                      const int& npol,
                                      double* emin_in,
                                      double* emax_in)
    : HamiltPW<T, Device>(pot_in, wfc_basis, p_kv, nlpp,ucell), ngk(p_kv->ngk)
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
    auto call_act = [&, this](const Operator<T, Device>* op, const bool& is_first_node) -> void {
        op->act(nbands, this->npwk_max, this->npol, psi_in, hpsi, this->ngk[op->get_ik()],  is_first_node);
    };

    ModuleBase::timer::tick("HamiltSdftPW", "hPsi");
    call_act(this->ops, true); // first node
    Operator<T, Device>* node((Operator<T, Device>*)this->ops->next_op);
    while (node != nullptr)
    {
        call_act(node, false); // other nodes
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

    hpsi_norm_op<Real, Device>()(this->ctx, nbands, npwk_max, npwk, Ebar, DeltaE, hpsi_norm, psi_in);
    ModuleBase::timer::tick("HamiltSdftPW", "hPsi_norm");
}

// template class HamiltSdftPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HamiltSdftPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
// template class HamiltSdftPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HamiltSdftPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt