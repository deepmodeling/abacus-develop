#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"

namespace hsolver {

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class HSolverPW: public HSolver<FPTYPE, Device>
{
  public:
    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in, wavefunc* pwf_in);

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    void solve(hamilt::Hamilt<FPTYPE, Device>* pHamilt,
               psi::Psi<std::complex<FPTYPE>, Device>& psi,
               elecstate::ElecState* pes,
               const std::string method_in,
               const bool skip_charge) override;

    virtual FPTYPE cal_hsolerror() override;
    virtual FPTYPE set_diagethr(const int istep, const int iter, const FPTYPE drho) override;
    virtual FPTYPE reset_diagethr(std::ofstream& ofs_running, const FPTYPE hsover_error, const FPTYPE drho) override;
  protected:
    void initDiagh(const psi::Psi<std::complex<FPTYPE>, Device>& psi_in);
    void endDiagh();
    void hamiltSolvePsiK(hamilt::Hamilt<FPTYPE, Device>* hm, psi::Psi<std::complex<FPTYPE>, Device>& psi, FPTYPE* eigenvalue);

    void updatePsiK(hamilt::Hamilt<FPTYPE, Device>* pHamilt,
                    psi::Psi<std::complex<FPTYPE>, Device>& psi,
                    const int ik);

    ModulePW::PW_Basis_K* wfc_basis = nullptr;
    wavefunc* pwf = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<FPTYPE> &h_diag, const int ik, const int npw);

    std::vector<FPTYPE> precondition;

    bool initialed_psi = false;

    Device * ctx = {};
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, psi::DEVICE_CPU>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, psi::DEVICE_CPU>;
    using castmem_2d_2h_op = psi::memory::cast_memory_op<double, FPTYPE, psi::DEVICE_CPU, psi::DEVICE_CPU>;
};

} // namespace hsolver

#endif