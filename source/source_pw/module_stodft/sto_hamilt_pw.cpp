#include "sto_hamilt_pw.h"
#include "source_base/timer.h"
#include "kernels/hpsi_norm_op.h"

template <typename T, typename Device>
StoHamiltPW<T, Device>::StoHamiltPW(elecstate::Potential* pot_in,
                                    ModulePW::PW_Basis_K* wfc_basis,
                                    K_Vectors* p_kv,
                                    pseudopot_cell_vnl* nlpp,
                                    const UnitCell* ucell,
                                    const int& npol,
                                    Real* emin_in,
                                    Real* emax_in)
    : hamilt::HamiltPW<T, Device>(pot_in, wfc_basis, p_kv, nlpp, nullptr, ucell, nullptr), ngk(p_kv->ngk)
{
    this->classname = "StoHamiltPW";
    this->npwk_max = wfc_basis->npwk_max;
    this->npol = npol;
    this->emin = emin_in;
    this->emax = emax_in;
}

template <typename T, typename Device>
void StoHamiltPW<T, Device>::hPsi(const T* psi_in, T* hpsi, const int& nbands)
{
    auto call_act = [&, this](const hamilt::Operator<T, Device>* op, const bool& is_first_node) -> void {
        op->act(nbands, this->npwk_max, this->npol, psi_in, hpsi, this->ngk[op->get_ik()],  is_first_node);
    };

    ModuleBase::timer::start("StoHamiltPW", "hPsi");
    call_act(this->ops, true); // first node
    hamilt::Operator<T, Device>* node((hamilt::Operator<T, Device>*)this->ops->next_op);
    while (node != nullptr)
    {
        call_act(node, false); // other nodes
        node = (hamilt::Operator<T, Device>*)(node->next_op);
    }
    ModuleBase::timer::end("StoHamiltPW", "hPsi");

    return;
}

template <typename T, typename Device>
void StoHamiltPW<T, Device>::hPsi_norm(const T* psi_in, T* hpsi_norm, const int& nbands)
{
    ModuleBase::timer::start("StoHamiltPW", "hPsi_norm");

    this->hPsi(psi_in, hpsi_norm, nbands);

    const int ik = this->ops->get_ik();
    const int npwk_max = this->npwk_max;
    const int npwk = this->ngk[ik];
    const Real emin = *this->emin;
    const Real emax = *this->emax;
    const Real Ebar = (emin + emax) / 2;
    const Real DeltaE = (emax - emin) / 2;

    hamilt::hpsi_norm_op<Real, Device>()(this->ctx, nbands, npwk_max, npwk, Ebar, DeltaE, hpsi_norm, psi_in);
    ModuleBase::timer::end("StoHamiltPW", "hPsi_norm");
}

template class StoHamiltPW<std::complex<float>, base_device::DEVICE_CPU>;
template class StoHamiltPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class StoHamiltPW<std::complex<float>, base_device::DEVICE_GPU>;
template class StoHamiltPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
