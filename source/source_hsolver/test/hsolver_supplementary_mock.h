#pragma once
#include "source_estate/elecstate_pw.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    // hamilt::MatrixBlock<double> temp{&(this->charge->rho[spin][0]), 1, this->charge->nrxx}; //
    // this->chr->get_nspin(), this->chr->get_nrxx()};
    return &(this->charge->rho[spin][0]);
}


void ElecState::init_nelec_spin()
{
    return;
}



void ElecState::init_ks(Charge* chg_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        const ModulePW::PW_Basis_Big* bigpw_in)
{
    return;
}

template <typename T, typename Device>
ElecStatePW<T, Device>::ElecStatePW(ModulePW::PW_Basis_K* wfc_basis_in,
                                    Charge* chg_in,
                                    K_Vectors* pkv_in,
                                    UnitCell* ucell_in,
                                    pseudopot_cell_vnl* ppcell_in,
                                    ModulePW::PW_Basis* rhopw_in,
                                    ModulePW::PW_Basis_Big* bigpw_in)
    : basis(wfc_basis_in)
{
}

template <typename T, typename Device>
ElecStatePW<T, Device>::~ElecStatePW()
{
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::psiToRho(const psi::Psi<T, Device>& psi)
{
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::cal_tau(const psi::Psi<T, Device>& psi)
{
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::cal_becsum(const psi::Psi<T, Device>& psi)
{
}

template <typename T, typename Device>
double ElecStatePW<T, Device>::get_spin_constrain_energy()
{
    return 0.0;
}

template class ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

Potential::~Potential()
{
}

void Potential::cal_v_eff(const Charge* const chg, const UnitCell* const ucell, ModuleBase::matrix& v_eff)
{
}

void Potential::cal_fixed_v(double* vl_pseudo)
{
}

} // namespace elecstate

#include "source_cell/klist.h"
