#ifndef NONLOCALPW_H
#define NONLOCALPW_H

#include "op_pw.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/nonlocal_workspace.h"

#include "source_pw/module_pwdft/vnl_pw.h"

namespace hamilt {

#ifndef NONLOCALTEMPLATE_H
#define NONLOCALTEMPLATE_H

template<class T> class Nonlocal : public T {};
// template<typename Real, typename Device = base_device::DEVICE_CPU>
// class Nonlocal : public OperatorPW<T, Device> {};

#endif

template<typename T, typename Device>
class Nonlocal<OperatorPW<T, Device>> : public OperatorPW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    Nonlocal(const int* isk_in,
             const pseudopot_cell_vnl* ppcell_in,
             const UnitCell* ucell_in,
             const ModulePW::PW_Basis_K* wfc_basis);

    template<typename T_in, typename Device_in = Device>
    explicit Nonlocal(const Nonlocal<OperatorPW<T_in, Device_in>>* nonlocal);

    virtual ~Nonlocal();

    virtual void init(const int ik_in)override;

    virtual void act(const int nbands,
        const int nbasis,
        const int npol,
        const T* tmpsi_in,
        T* tmhpsi,
        const int ngk_ik = 0,
        const bool is_first_node = false)const override;

    const int *get_isk() const {return this->isk;}
    const pseudopot_cell_vnl *get_ppcell() const {return this->ppcell;}
    const UnitCell *get_ucell() const {return this->ucell;}
    /** @brief Return the borrowed plane-wave basis and pool communicator. */
    const ModulePW::PW_Basis_K* get_wfcpw() const { return this->wfcpw; }
    T* get_vkb() const
    {
        return this->vkb;
    }
    T* get_becp() const
    {
        return this->workspace_.get_becp();
    }

  private:
    mutable NonlocalWorkspace<T, Device> workspace_;

    const int* isk = nullptr;

    const pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell* ucell = nullptr;

    const ModulePW::PW_Basis_K* wfcpw = nullptr;

    mutable T *vkb = nullptr;
    Device* ctx = {};

};

} // namespace hamilt

#endif
