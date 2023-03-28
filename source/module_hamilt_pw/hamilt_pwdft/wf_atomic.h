#ifndef WF_ATOMIC_H
#define WF_ATOMIC_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/complexmatrix.h"
#include "wf_igk.h"
#include "module_psi/psi.h"
#include "module_basis/module_pw/pw_basis_k.h"

class WF_atomic : public WF_igk
{
	public:

    WF_atomic();
    ~WF_atomic();

	ModuleBase::realArray table_local;//mohan add 2009-09-10

    ModuleBase::ComplexMatrix *evc = nullptr;  // wavefunctions in the PW basis
    //temporary psi for new code
    psi::Psi<std::complex<double>>* psi = nullptr;
    void evc_transform_psi();
    void psi_transform_evc();

    ModuleBase::ComplexMatrix *wanf2 = nullptr; // wannier functions in the PW basis
    
    int pw_seed; //random seed for wave functions qianrui add 2021-8-13

    void init_at_1(void);// from init_at_1.f90

    void print_PAOs(void)const;

    public: //template change to public, will be refactor later. added by zhengdy 20230302
    int *irindex = nullptr;

    void atomic_wfc
    (
        const int ik,
        const int np,
        const int lmax_wfc,
        ModuleBase::ComplexMatrix &wfcatom,
        const ModuleBase::realArray &table_q,
        const int &table_dimension,
        const double &dq
    )const;

    //==================================
    // Calculate random wave functions
    // as trial wave functions
    //==================================
    void atomicrandom(ModuleBase::ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis)const;

    void random(std::complex<double> *psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis);
    void random(std::complex<float> *psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis);

    template<typename FPTYPE>
    void random_t(std::complex<FPTYPE> *psi,const int iw_start,const int iw_end,const int ik, ModulePW::PW_Basis_K* wfc_basis);

    void check_evc()const;
#ifdef __MPI
	void stick_to_pool(double *stick, const int &ir, double *out, ModulePW::PW_Basis_K* wfc_basis) const;
    void stick_to_pool(float *stick, const int &ir, float *out, ModulePW::PW_Basis_K* wfc_basis) const;
#endif

};

#endif 
