#ifndef OP_PW_EKIN_TD_H
#define OP_PW_EKIN_TD_H

#include "op_pw.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_basis/module_pw/pw_basis_k.h"

namespace hamilt
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class TDEkineticPW : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    /** @brief Bind the vector potential in Hartree atomic units before init(). */
    void set_A_ha(const ModuleBase::Vector3<double>& A) { A_ha_ = A; }
    TDEkineticPW(const ModulePW::PW_Basis_K* wfc_basis_in);

    virtual ~TDEkineticPW();

    void init(const int ik) override;

    virtual void act(const int nbands,
                     const int nbasis,
                     const int npol,
                     const T* tmpsi_in,
                     T* tmhpsi,
                     const int ngk_ik = 0,
                     const bool is_first_node = false) const override;

  private:
    const ModulePW::PW_Basis_K* wfc_basis = nullptr;
    ModuleBase::Vector3<double> A_ha_;
    Device* ctx = {};
    ct::Tensor factor_;
};

} // namespace hamilt

#endif // OP_PW_EKIN_TD_H
