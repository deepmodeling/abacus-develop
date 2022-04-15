#ifndef WF_ATOMIC_H
#define WF_ATOMIC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/complexmatrix.h"
#include "wf_igk.h"

class WF_atomic : public WF_igk
{
	public:

    WF_atomic();
    ~WF_atomic();

	ModuleBase::realArray table_local;//mohan add 2009-09-10

    ModuleBase::ComplexMatrix *evc;  // wavefunctions in the PW basis

    ModuleBase::ComplexMatrix *wanf2; // wannier functions in the PW basis
    
    int pw_seed; //random seed for wave functions qianrui add 2021-8-13

    void init_at_1(void);// from init_at_1.f90

    void print_PAOs(void)const;

	protected:

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
    void random(ModuleBase::ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const;
    void atomicrandom(ModuleBase::ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const;

    void check_psi(const ModuleBase::ComplexMatrix *psi)const;

};
#endif 
