#ifndef __VEFFPW
#define __VEFFPW

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt
{

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template<class T> class Veff : public T 
{};

#endif

template<>
class Veff<OperatorPW> : public OperatorPW
{
    public:
    Veff(
        const int* isk_in,
        const ModuleBase::matrix* veff_in,
        ModulePW::PW_Basis_K* wfcpw_in
    );

    virtual ~Veff(){};

    virtual void act
    (
        const psi::Psi<std::complex<double>> *psi_in, 
        const int n_npwx, 
        const std::complex<double>* tmpsi_in, 
        std::complex<double>* tmhpsi
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