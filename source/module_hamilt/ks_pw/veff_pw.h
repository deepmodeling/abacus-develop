#ifndef VEFFPW_H
#define VEFFPW_H

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt {

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template<class T> class Veff : public T
{};
// template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
// class Veff : public OperatorPW<FPTYPE, Device> {};

#endif

template<typename FPTYPE, typename Device>
class Veff<OperatorPW<FPTYPE, Device>> : public OperatorPW<FPTYPE, Device>
{
  public:
    Veff(const int* isk_in,const ModuleBase::matrix* veff_in,ModulePW::PW_Basis_K* wfcpw_in);

    virtual ~Veff(){};

    virtual void act (
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi
    )const override;

  private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    const int* isk = nullptr;

    const ModuleBase::matrix* veff = nullptr;

    ModulePW::PW_Basis_K* wfcpw = nullptr;
};

} // namespace hamilt

#endif