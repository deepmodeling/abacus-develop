#ifndef EKINETICPW_H
#define EKINETICPW_H

#include "operator_pw.h"
#include "module_hamilt/include/ekinetic.h"

namespace hamilt {

// Not needed anymore
#ifndef __EKINETICTEMPLATE
#define __EKINETICTEMPLATE

// template<class T> class Ekinetic : public T {};
template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
class Ekinetic : public OperatorPW<FPTYPE, Device> {};

#endif

// template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
// class Ekinetic : public OperatorPW<FPTYPE, Device>
template<typename FPTYPE, typename Device>
class Ekinetic<OperatorPW<FPTYPE, Device>> : public OperatorPW<FPTYPE, Device>
{
  public:
    Ekinetic(
        FPTYPE tpiba2_in, 
        const FPTYPE* gk2_in,
        const int gk2_row, 
        const int gk2_col);

    virtual ~Ekinetic();

    virtual void act(
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi)const override;

  private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    FPTYPE tpiba2 = 0.0;

#if ((defined __CUDA) || (defined __ROCM))
    FPTYPE* gk2 = nullptr;
#else
    const FPTYPE* gk2 = nullptr;
#endif
    int gk2_row = 0;
    int gk2_col = 0;

    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};

    using ekinetic_op = ekinetic_pw_op<FPTYPE, Device>;
    using resize_memory_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delete_memory_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using synchronize_memory_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
};

} // namespace hamilt

#endif