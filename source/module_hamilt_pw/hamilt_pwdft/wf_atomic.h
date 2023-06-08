#ifndef WF_ATOMIC_H
#define WF_ATOMIC_H

#include "module_base/complexmatrix.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/realarray.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_psi/psi.h"
#include "structure_factor.h"

class WF_atomic
{
	public:

    WF_atomic();
    ~WF_atomic();
    int npwx;
    int npw;
    // ModuleBase::IntArray igk;
#ifdef __CUDA
    double *d_g2kin;
#endif

    ModuleBase::realArray table_local;//mohan add 2009-09-10

    //temporary psi for new code
    psi::Psi<std::complex<double>>* psi = nullptr;

    ModuleBase::ComplexMatrix *wanf2 = nullptr; // wannier functions in the PW basis

    void init_at_1(Structure_Factor *sf_in); // from init_at_1.f90

    void print_PAOs(void)const;

    public: //template change to public, will be refactor later. added by zhengdy 20230302
    int *irindex = nullptr;

    void atomic_wfc(const int ik,
                    const int np,
                    const int lmax_wfc,
                    const ModulePW::PW_Basis_K* wfc_basis,
                    ModuleBase::ComplexMatrix& wfcatom,
                    const ModuleBase::realArray& table_q,
                    const int& table_dimension,
                    const double& dq) const;

    //==================================
    // Calculate random wave functions
    // as trial wave functions
    //==================================
    void atomicrandom(ModuleBase::ComplexMatrix& psi,
                      const int iw_start,
                      const int iw_end,
                      const int ik,
                      const ModulePW::PW_Basis_K* wfc_basis) const;

    void random(std::complex<double>* psi,
                const int iw_start,
                const int iw_end,
                const int ik,
                const ModulePW::PW_Basis_K* wfc_basis);
    void random(std::complex<float>* psi,
                const int iw_start,
                const int iw_end,
                const int ik,
                const ModulePW::PW_Basis_K* wfc_basis);

    template <typename FPTYPE>
    void random_t(std::complex<FPTYPE>* psi,
                  const int iw_start,
                  const int iw_end,
                  const int ik,
                  const ModulePW::PW_Basis_K* wfc_basis);

#ifdef __MPI
    void stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const;
    void stick_to_pool(float* stick, const int& ir, float* out, const ModulePW::PW_Basis_K* wfc_basis) const;
#endif
  private:
    Structure_Factor *psf;
};

#endif 
