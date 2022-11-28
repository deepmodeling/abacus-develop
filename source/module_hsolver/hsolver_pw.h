#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_pw/pw_basis_k.h"

namespace hsolver {

template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class HSolverPW: public HSolver<FPTYPE, Device>
{
  public:
    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in);

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    void solve(hamilt::Hamilt<FPTYPE, Device>* pHamilt, psi::Psi<std::complex<FPTYPE>, Device>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    virtual FPTYPE cal_hsolerror() override;
    virtual FPTYPE set_diagethr(const int istep, const int iter, const FPTYPE drho) override;
    virtual FPTYPE reset_diagethr(std::ofstream& ofs_running, const FPTYPE hsover_error, const FPTYPE drho) override;
  protected:
    void initDiagh();
    void endDiagh();
    void hamiltSolvePsiK(hamilt::Hamilt<FPTYPE, Device>* hm, psi::Psi<std::complex<FPTYPE>, Device>& psi, FPTYPE* eigenvalue);

    void updatePsiK(hamilt::Hamilt<FPTYPE, Device>* pHamilt, psi::Psi<std::complex<FPTYPE>, Device>& psi, const int ik);

    ModulePW::PW_Basis_K* wfc_basis = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<FPTYPE> &h_diag, const int ik, const int npw);

    std::vector<FPTYPE> precondition;

    bool initialed_psi = false;

    Device *ctx = {};
};

} // namespace hsolver

#endif