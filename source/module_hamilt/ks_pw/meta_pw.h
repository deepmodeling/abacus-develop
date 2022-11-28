#ifndef METAPW_H
#define METAPW_H

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"

namespace hamilt {

#ifndef __METATEMPLATE
#define __METATEMPLATE

template<class T> class Meta : public T
{};
// template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
// class Meta : public OperatorPW<FPTYPE, Device> {};

#endif

template<typename FPTYPE, typename Device>
class Meta<OperatorPW<FPTYPE, Device>> : public OperatorPW<FPTYPE, Device>
{
    public:
    Meta(FPTYPE tpiba2_in, const int* isk_in, const ModuleBase::matrix* vk, ModulePW::PW_Basis_K* wfcpw);

    template<typename T_in, typename Device_in = Device>
    explicit Meta(const Meta<OperatorPW<T_in, Device_in>>* meta);

    virtual ~Meta(){};

    virtual void act(
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi
    )const override;

    // denghui added for copy constructor at 20221105
    FPTYPE get_tpiba() const {return this->tpiba;}
    const int * get_isk() const {return this->isk;}
    const ModuleBase::matrix* get_vk() const {return this->vk;}
    ModulePW::PW_Basis_K* get_wfcpw() const {return this->wfcpw;}

    private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    FPTYPE tpiba = 0.0;

    const int* isk = nullptr;

    const ModuleBase::matrix* vk = nullptr;

    ModulePW::PW_Basis_K* wfcpw = nullptr;
};

} // namespace hamilt

#endif