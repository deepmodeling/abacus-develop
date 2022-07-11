#ifndef HSOLVERPW_H
#define HSOLVERPW_H

#include "hsolver.h"
#include "module_pw/pw_basis_k.h"

namespace hsolver
{

class HSolverPW : public HSolver
{
  public:
    HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in);

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/

    void solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge) override;

    virtual double cal_hsolerror() override;
    virtual double set_diagethr(const int istep, const int iter, const double drho) override;
    virtual double reset_diagethr(std::ofstream& ofs_running, const double hsover_error, const double drho) override;
  protected:
    void initDiagh();
    void endDiagh();
    void hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue);

    void updatePsiK(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, const int ik);

    ModulePW::PW_Basis_K* wfc_basis = nullptr;

    // calculate the precondition array for diagonalization in PW base
    void update_precondition(std::vector<double> &h_diag, const int ik, const int npw);

    std::vector<double> precondition;

    bool initialed_psi = false;
};

} // namespace hsolver

#endif