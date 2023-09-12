#ifndef EKINETICPW_H
#define EKINETICPW_H

#include "operator_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/ekinetic_op.h"

namespace hamilt {

// Not needed anymore
#ifndef __EKINETICTEMPLATE
#define __EKINETICTEMPLATE

template<class T> class Ekinetic : public T
{};
// template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
// class Ekinetic : public OperatorPW<FPTYPE, Device> {};

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

    template<typename T_in, typename Device_in = Device>
    explicit Ekinetic(const Ekinetic<OperatorPW<T_in, Device_in>>* ekinetic);

    virtual ~Ekinetic();

    virtual void act(const int nbands,
        const int nbasis,
        const int npol,
        const std::complex<FPTYPE>* tmpsi_in,
        std::complex<FPTYPE>* tmhpsi,
        const int ngk_ik = 0)const override;

    // denghuilu added for copy construct at 20221105
    int get_gk2_row() const {return this->gk2_row;}
    int get_gk2_col() const {return this->gk2_col;}
    FPTYPE get_tpiba2() const {return this->tpiba2;}
    const FPTYPE* get_gk2() const {return this->gk2;}
    Device* get_ctx() const {return this->ctx;}

  private:

    FPTYPE tpiba2 = 0.0;
    const FPTYPE* gk2 = nullptr;
    int gk2_row = 0;
    int gk2_col = 0;

    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    psi::AbacusDevice_t device = {};

    using ekinetic_op = ekinetic_pw_op<FPTYPE, Device>;
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using syncmem_var_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
};

} // namespace hamilt

#endif