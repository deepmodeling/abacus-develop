#ifndef VEFFPW_H
#define VEFFPW_H

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/veff_op.h"

#include <module_base/macros.h>

namespace hamilt {

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template<class T> class Veff : public T {};
// template<typename Real, typename Device = psi::DEVICE_CPU>
// class Veff : public OperatorPW<T, Device> {};

#endif

template<typename T, typename Device>
class Veff<OperatorPW<T, Device>> : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    Veff(const int* isk_in,
         const Real* veff_in,
         const int veff_row,
         const int veff_col,
         const ModulePW::PW_Basis_K* wfcpw_in);

    template<typename T_in, typename Device_in = Device>
    explicit Veff(const Veff<OperatorPW<T_in, Device_in>>* veff);

    virtual ~Veff();

    virtual void act(const int nbands,
        const int nbasis,
        const int npol,
        const T* tmpsi_in,
        T* tmhpsi,
        const int ngk_ik = 0)const override;

    // denghui added for copy constructor at 20221105
    const Real *get_veff() const {return this->veff;}
    int get_veff_col() const {return this->veff_col;}
    int get_veff_row() const { return this->veff_row; }
    const int *get_isk() const {return isk;}
    const ModulePW::PW_Basis_K* get_wfcpw() const
    {
        return this->wfcpw;
    }

  private:

    const int* isk = nullptr;

    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};

    int veff_col = 0;
    int veff_row = 0;
    const Real *veff = nullptr, *h_veff = nullptr, *d_veff = nullptr;
    T *porter = nullptr;
    T *porter1 = nullptr;
    psi::AbacusDevice_t device = {};
    using veff_op = veff_pw_op<Real, Device>;


    using resmem_complex_op = psi::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<T, Device>;
};

} // namespace hamilt

#endif