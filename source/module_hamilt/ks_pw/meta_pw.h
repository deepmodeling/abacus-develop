#ifndef __METAPW
#define __METAPW

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt
{

#ifndef __METATEMPLATE
#define __METATEMPLATE

template<class T> class Meta : public T
{};

#endif

template<>
class Meta<OperatorPW> : public OperatorPW
{
    public:
    Meta(
        double tpiba2_in,
        const int* isk_in,
        const ModuleBase::matrix* vk,
        ModulePW::PW_Basis_K* wfcpw
    );

    virtual ~Meta(){};

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

    double tpiba = 0.0;

    const int* isk = nullptr;

    const ModuleBase::matrix* vk = nullptr;

    ModulePW::PW_Basis_K* wfcpw = nullptr;
};

} // namespace hamilt

#endif