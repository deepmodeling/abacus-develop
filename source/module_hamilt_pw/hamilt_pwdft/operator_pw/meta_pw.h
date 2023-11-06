#ifndef METAPW_H
#define METAPW_H

#include "operator_pw.h"
#include "module_base/matrix.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hamilt_pw/hamilt_pwdft/kernels/meta_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"

#include <module_base/macros.h>

namespace hamilt {

#ifndef __METATEMPLATE
#define __METATEMPLATE

template<class T> class Meta : public T {};
// template<typename Real, typename Device = psi::DEVICE_CPU>
// class Meta : public OperatorPW<T, Device> {};

#endif

template<typename T, typename Device>
class Meta<OperatorPW<T, Device>> : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;
    public:
      Meta(Real tpiba2_in,
           const int* isk_in,
           const Real* vk_in,
           const int vk_row,
           const int vk_col,
           const ModulePW::PW_Basis_K* wfcpw);

      template <typename T_in, typename Device_in = Device>
      explicit Meta(const Meta<OperatorPW<T_in, Device_in>>* meta);

      virtual ~Meta();

      virtual void act(const int nbands,
          const int nbasis,
          const int npol,
          const T* tmpsi_in,
          T* tmhpsi,
          const int ngk = 0)const override;

      // denghui added for copy constructor at 20221105
      Real get_tpiba() const
      {
          return this->tpiba;
      }
    const int * get_isk() const {return this->isk;}
    const Real* get_vk() const {return this->vk;}
    int get_vk_row() const {return this->vk_row;}
    int get_vk_col() const {return this->vk_col;}
    const ModulePW::PW_Basis_K* get_wfcpw() const
    {
          return this->wfcpw;
    }

    private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    mutable int vk_row = 0;
    mutable int vk_col = 0;

    Real tpiba = 0.0;

    const int* isk = nullptr;

    const Real * vk = nullptr;

    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    Device* ctx = {};
    psi::DEVICE_CPU* cpu_ctx = {};
    T *porter = nullptr;
    using meta_op = meta_pw_op<Real, Device>;
    using vector_mul_vector_op = hsolver::vector_mul_vector_op<T, Device>;
    using resmem_complex_op = psi::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = psi::memory::delete_memory_op<T, Device>;
};

} // namespace hamilt

#endif