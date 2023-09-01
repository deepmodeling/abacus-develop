#ifndef OPERATORPW_H
#define OPERATORPW_H
#include"module_hamilt_general/operator.h"

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

    std::string classname = "";
    using syncmem_complex_op = psi::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, Device>;
};

}//end namespace hamilt

#endif