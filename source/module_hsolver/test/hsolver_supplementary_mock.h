#pragma once
#include "module_elecstate/elecstate_pw.h"

namespace elecstate
{

const double* ElecState::getRho(int spin) const
{
    // hamilt::MatrixBlock<double> temp{&(this->charge->rho[spin][0]), 1, this->charge->nrxx}; //
    // this->chr->get_nspin(), this->chr->get_nrxx()};
    return &(this->charge->rho[spin][0]);
}

void ElecState::fixed_weights(const std::vector<double>& ocp_kb, const int& nbands, const double& nelec)
{
    return;
}

void ElecState::init_nelec_spin()
{
    return;
}

void ElecState::calculate_weights()
{
    return;
}

void ElecState::calEBand()
{
    return;
}

void ElecState::print_band(const int& ik, const int& printe, const int& iter)
{
    return;
}

void ElecState::print_eigenvalue(std::ofstream& ofs)
{
    return;
}

void ElecState::init_scf(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Grid& pgrid,
                         const ModuleBase::ComplexMatrix& strucfac,
                         const bool*,
                         ModuleSymmetry::Symmetry&,
                         const void*)
{
    return;
}

void ElecState::init_ks(Charge* chg_in, // pointer for class Charge
                        const K_Vectors* klist_in,
                        int nk_in,
                        ModulePW::PW_Basis* rhopw_in,
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
                                    ModulePW::PW_Basis* rhodpw_in,
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

// mock of Stochastic_WF
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
template <typename T, typename Device>
Stochastic_WF<T, Device>::Stochastic_WF()
{
    chiortho = nullptr;
    chi0 = nullptr;
    shchi = nullptr;
    nchip = nullptr;
}

template <typename T, typename Device>
Stochastic_WF<T, Device>::~Stochastic_WF()
{
    delete[] chi0;
    delete[] shchi;
    delete[] chiortho;
    delete[] nchip;
}

template <typename T, typename Device>
void Stochastic_WF<T, Device>::init(K_Vectors* p_kv, const int npwx_in)
{
    /*chi0 = new ModuleBase::ComplexMatrix[nks_in];
    shchi = new ModuleBase::ComplexMatrix[nks_in];
    chiortho = new ModuleBase::ComplexMatrix[nks_in];
    nchip = new int[nks_in];
    this->nks = nks_in;*/
}

#include "module_cell/klist.h"
K_Vectors::K_Vectors()
{
}
K_Vectors::~K_Vectors()
{
}
