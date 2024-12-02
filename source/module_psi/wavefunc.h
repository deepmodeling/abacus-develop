#ifndef WAVEFUNC_H
#define WAVEFUNC_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_hamilt_general/hamilt.h"
#include "module_psi/wf_atomic.h"

class wavefunc : public WF_atomic
{
  public:
    wavefunc();
    ~wavefunc();

    // allocate memory
    psi::Psi<std::complex<double>>* allocate(const int nkstot, const int nks, const int* ngk, const int npwx);

    int nkstot = 0; // total number of k-points for all pools

    void wfcinit(psi::Psi<std::complex<double>>* psi_in, ModulePW::PW_Basis_K* wfc_basis);

    int get_starting_nw(void) const;

    void init_after_vc(const int nks); // LiuXh 20180515
};

namespace hamilt
{

void diago_PAO_in_pw_k2(const int& ik,
                        psi::Psi<std::complex<float>>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<float>>* phm_in = nullptr);
void diago_PAO_in_pw_k2(const int& ik,
                        psi::Psi<std::complex<double>>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<double>>* phm_in = nullptr);
void diago_PAO_in_pw_k2(const int& ik,
                        ModuleBase::ComplexMatrix& wvf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        wavefunc* p_wf);

template <typename FPTYPE, typename Device>
void diago_PAO_in_pw_k2(const Device* ctx,
                        const int& ik,
                        psi::Psi<std::complex<FPTYPE>, Device>& wvf,
                        ModulePW::PW_Basis_K* wfc_basis,
                        wavefunc* p_wf,
                        const ModuleBase::realArray& tab_at,
                        const int& lmaxkb,
                        hamilt::Hamilt<std::complex<FPTYPE>, Device>* phm_in = nullptr);
} // namespace hamilt

#endif // wavefunc
