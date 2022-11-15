#ifndef VEFFPW_H
#define VEFFPW_H

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis_k.h"
#include "module_hamilt/include/veff.h"

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

    template<typename T_in, typename Device_in = Device>
    explicit Veff(const Veff<OperatorPW<T_in, Device_in>>* veff);

    virtual ~Veff();

    virtual void act (
        const psi::Psi<std::complex<FPTYPE>, Device> *psi_in, 
        const int n_npwx, 
        const std::complex<FPTYPE>* tmpsi_in, 
        std::complex<FPTYPE>* tmhpsi
    )const override;

    // denghui added for copy constructor at 20221105
    FPTYPE *get_veff() const {return this->veff;}
    int get_veff_col() const {return this->veff_col;}
    int get_veff_row() const {return this->veff_row;}
    int get_npol() const {return this->npol;}
    const int *get_isk() const {return isk;}
    ModulePW::PW_Basis_K* get_wfcpw() const {return this->wfcpw;}

  private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    const int* isk = nullptr;


    ModulePW::PW_Basis_K* wfcpw = nullptr;

    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};

    int veff_col = 0;
    int veff_row = 0;
    FPTYPE *veff = nullptr;
    FPTYPE *d_veff = nullptr;
    std::complex<FPTYPE> *porter = nullptr;
    std::complex<FPTYPE> *porter1 = nullptr;
    psi::AbacusDevice_t device = {};
    using veff_op = veff_pw_op<FPTYPE, Device>;

    using resize_memory_double_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delete_memory_double_op = psi::memory::delete_memory_op<FPTYPE, Device>;

    using resize_memory_complex_op = psi::memory::resize_memory_op<std::complex<FPTYPE>, Device>;
    using delete_memory_complex_op = psi::memory::delete_memory_op<std::complex<FPTYPE>, Device>;
    using syncmem_double_h2d_op = psi::memory::synchronize_memory_op<FPTYPE, Device, psi::DEVICE_CPU>;
};

} // namespace hamilt

#endif