#ifndef OP_PW_NL_TD_H
#define OP_PW_NL_TD_H

#include "op_pw.h"
#include "source_pw/module_pwdft/nonlocal_workspace.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace hamilt
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class TDNonlocalPW : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    /** @brief Bind the vector potential in Hartree atomic units before init(). */
    void set_A_ha(const ModuleBase::Vector3<double>& A) { A_ha_ = A; }
    TDNonlocalPW(const int* isk_in, const pseudopot_cell_vnl* ppcell_in, const UnitCell* ucell_in, const ModulePW::PW_Basis_K* wfc_basis);

    virtual ~TDNonlocalPW();
    virtual void init(const int ik_in) override;

    virtual void act(const int nbands,
                     const int nbasis,
                     const int npol,
                     const T* tmpsi_in,
                     T* tmhpsi,
                     const int ngk_ik = 0,
                     const bool is_first_node = false) const override;

  private:
    mutable NonlocalWorkspace<T, Device> workspace_;

    const int* isk = nullptr;
    const pseudopot_cell_vnl* ppcell = nullptr;
    const UnitCell* ucell = nullptr;
    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    mutable T* vkb_td = nullptr; // Time-dependent projector cache.

    ModuleBase::Vector3<double> A_ha_;
    Device* ctx = {};
    using resmem_complex_op = base_device::memory::resize_memory_op<T, Device>;
    using delmem_complex_op = base_device::memory::delete_memory_op<T, Device>;

};

} // namespace hamilt
#endif
