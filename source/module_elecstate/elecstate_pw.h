#ifndef ELECSTATEPW_H
#define ELECSTATEPW_H

#include "elecstate.h"
#include "module_pw/pw_basis_k.h"
#include "module_elecstate/include/elecstate_multi_device.h"

namespace elecstate
{

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class ElecStatePW : public ElecState
{
  public:
    ElecStatePW(
      ModulePW::PW_Basis_K *wfc_basis_in, 
      Charge* chg_in, 
      K_Vectors *pkv_in
    );
    // void init(Charge* chg_in):charge(chg_in){} override;

    ~ElecStatePW();
    // interface for HSolver to calculate rho from Psi
    virtual void psiToRho(const psi::Psi<std::complex<FPTYPE>, Device>& psi);
    // return current electronic density rho, as a input for constructing Hamiltonian
    // const double* getRho(int spin) const override;

    // update charge density for next scf step
    // void getNewRho() override;

  protected:
    ModulePW::PW_Basis_K *basis;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    void updateRhoK(const psi::Psi<std::complex<FPTYPE>, Device>& psi); // override;
    // sum over all pools for rho and ebands
    void parallelK();
    // calcualte rho for each k
    void rhoBandK(const psi::Psi<std::complex<FPTYPE>, Device>& psi);

    void init_rho_data();

    Device * ctx = {};
    psi::DEVICE_CPU * cpu_ctx = {};
    bool init_rho = false;
    FPTYPE ** rho = nullptr, ** kin_r = nullptr;
    FPTYPE * rho_data = nullptr, * kin_r_data = nullptr;
    std::complex<FPTYPE> *wfcr = nullptr, *wfcr_another_spin = nullptr;

    using elecstate_pw_op = elecstate::elecstate_pw_op<FPTYPE, Device>;

    using setmem_var_op = psi::memory::set_memory_op<FPTYPE, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_d2h_op = psi::memory::synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, Device>;

    using resmem_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
};

} // namespace elecstate

#endif