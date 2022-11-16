#ifndef OPERATORPW_H
#define OPERATORPW_H
#include"module_hamilt/operator.h"

namespace hamilt {
template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class OperatorPW : public Operator<std::complex<FPTYPE>, Device>
{
    public:
    virtual ~OperatorPW();
    
    //in PW code, different operators donate hPsi independently
    //run this->act function for the first operator and run all act() for other nodes in chain table 
    using hpsi_info = typename hamilt::Operator<std::complex<FPTYPE>, Device>::hpsi_info;
    virtual hpsi_info hPsi(hpsi_info& input)const;
    //main function which should be modified in Operator for PW base
    //do operation : |hpsi_choosed> = V|psi_choosed>
    //V is the target operator act on choosed psi, the consequence should be added to choosed hpsi
    virtual void act(
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi)const;

    std::string classname = "";
};

}//end namespace hamilt

#endif