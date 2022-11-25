//==========================================================
// AUTHOR : wangjp
// Data :2009-04
// Last Update:
//
// 09-05-10 modify SchmitOrth() diag_zhegvx() as static
// member function
//==========================================================

#ifndef DIAGODAVID_H
#define DIAGODAVID_H

#include "diagh.h"
#include "module_base/complexmatrix.h"
#include "src_pw/structure_factor.h"

namespace hsolver {

template<typename FPTYPE = double, typename Device = psi::DEVICE_CPU>
class DiagoDavid : public DiagH<FPTYPE, Device>
{
  public:
    DiagoDavid(const FPTYPE* precondition_in);

    // this is the override function diag() for CG method
    void diag(hamilt::Hamilt<FPTYPE, Device>* phm_in, psi::Psi<std::complex<FPTYPE>, Device>& phi, FPTYPE* eigenvalue_in);

    static int PW_DIAG_NDIM;

  private:
    int test_david = 0;

    /// record for how many bands not have convergence eigenvalues
    int notconv = 0;

    /// device type of psi
    psi::AbacusDevice_t device = {};
    psi::DEVICE_CPU * ctx = {};

    void cal_grad(hamilt::Hamilt<FPTYPE, Device>* phm_in,
                  const int& npw,
                  const int& nbase,
                  const int& notconv,
                  psi::Psi<std::complex<FPTYPE>, Device>& basis,
                  ModuleBase::ComplexMatrix& hp,
                  ModuleBase::ComplexMatrix& sp,
                  const ModuleBase::ComplexMatrix& vc,
                  const int* unconv,
                  const FPTYPE* en);

    void cal_elem(const int& npw,
                  int& nbase,
                  const int& notconv,
                  const psi::Psi<std::complex<FPTYPE>, Device>& basis,
                  const ModuleBase::ComplexMatrix& hp,
                  const ModuleBase::ComplexMatrix& sp,
                  ModuleBase::ComplexMatrix& hc,
                  ModuleBase::ComplexMatrix& sc);

    void refresh(const int& npw,
                 const int& nband,
                 int& nbase,
                 const FPTYPE* en,
                 const psi::Psi<std::complex<FPTYPE>, Device>& psi,
                 psi::Psi<std::complex<FPTYPE>, Device>& basis,
                 ModuleBase::ComplexMatrix& hp,
                 ModuleBase::ComplexMatrix& sp,
                 ModuleBase::ComplexMatrix& hc,
                 ModuleBase::ComplexMatrix& sc,
                 ModuleBase::ComplexMatrix& vc);

    void cal_err(const int& npw,
                 const int& nband,
                 const int& nbase,
                 const ModuleBase::ComplexMatrix& vc,
                 const ModuleBase::ComplexMatrix& hp,
                 const psi::Psi<std::complex<FPTYPE>, Device>& basis,
                 const FPTYPE* en,
                 std::complex<FPTYPE>* respsi);

    void SchmitOrth(const int& npw,
                    const int n_band,
                    const int m,
                    psi::Psi<std::complex<FPTYPE>, Device>& psi,
                    const ModuleBase::ComplexMatrix& spsi,
                    std::complex<FPTYPE>* lagrange_m,
                    const int mm_size,
                    const int mv_size);
    void planSchmitOrth(
                    const int nband,
                    int* pre_matrix_mm_m,
                    int* pre_matrix_mv_m);

    void diag_zhegvx(const int& n,
                     const int& m,
                     const ModuleBase::ComplexMatrix& hc,
                     const ModuleBase::ComplexMatrix& sc,
                     const int& ldh,
                     FPTYPE* e,
                     ModuleBase::ComplexMatrix& vc);

    void diag_mock(hamilt::Hamilt<FPTYPE, Device>* phm_in, psi::Psi<std::complex<FPTYPE>, Device>& psi, FPTYPE* eigenvalue_in);

    const FPTYPE* precondition = nullptr;

    using hpsi_info = typename hamilt::Operator<std::complex<FPTYPE>, Device>::hpsi_info;
};

} // namespace hsolver

#endif
