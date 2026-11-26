#include "socket_driver.h"

#include "source_base/parallel_grid.h"
#include "source_cell/stru_fac.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver_ks_pw.h"
#include "source_estate/electronic_control.h"
#include "source_hamilt/module_xc/xc_functional.h"

#ifdef __LCAO
#include "source_esolver/esolver_ks_lcao.h"
#ifdef __EXX
#include "source_lcao/module_ri/exx_lri_interface.h"
#endif
#endif

namespace
{
// Reuse the native solver lifecycle; only a bound electronic request changes
// convergence handling or collects feedback after the native nuclear forces.
template <typename Solver>
class SocketKS : public Solver, public SocketCP::ElectronicSlot
{
  public:
    void cal_force(BaseCell &basecell, ModuleBase::matrix &force) override
    {
        Solver::cal_force(basecell, force);
        if (this->electronic_step != nullptr)
        {
            basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
            const UnitCell &cell = static_cast<const UnitCell &>(basecell);
            elecstate::complete_electronic_step(
                cell, this->pelec, this->pw_rhod, this->solvent, *this->inp_, this->conv_esolver,
                *this->electronic_step);
        }
    }

  protected:
    void iter_finish(
        UnitCell &cell,
        const int istep,
        int &iter,
        bool &converged) override
    {
        Solver::iter_finish(cell, istep, iter, converged);
        if (this->electronic_step != nullptr && this->oscillate_esolver)
        {
            converged = false;
        }
    }
};

class SocketPW : public SocketKS<ModuleESolver::ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>>
{
  protected:
    void before_scf(UnitCell &cell, const int istep) override
    {
        const bool active = this->electronic_step != nullptr;
        if (active)
        {
            this->conv_esolver = false;
            const int functional_type = XC_Functional::get_func_type();
            validate_electronic_step_exx(functional_type, this->general_exx_info_.cal_exx);
            elecstate::prepare_electronic_step(cell, this->pelec, *this->electronic_step, *this->inp_);
        }
        this->exx_helper->begin_electronic_step(active);
        ModuleESolver::ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>::before_scf(cell, istep);
    }
};

#ifdef __LCAO
template <typename TK>
class SocketLCAO : public SocketKS<ModuleESolver::ESolver_KS_LCAO<TK, double>>
{
  protected:
    void before_scf(UnitCell &cell, const int istep) override
    {
        const bool active = this->electronic_step != nullptr;
        if (active)
        {
            this->conv_esolver = false;
            elecstate::prepare_electronic_step(cell, this->pelec, *this->electronic_step, *this->inp_);
        }
#ifdef __EXX
        if (this->exx_nao.exd)
        {
            this->exx_nao.exd->begin_electronic_step(active);
        }
        if (this->exx_nao.exc)
        {
            this->exx_nao.exc->begin_electronic_step(active);
        }
#endif
        ModuleESolver::ESolver_KS_LCAO<TK, double>::before_scf(cell, istep);
    }
};
#endif
} // namespace

ModuleESolver::ESolver *SocketCP::create_pw_socket_solver() { return new SocketPW(); }

#ifdef __LCAO
ModuleESolver::ESolver *SocketCP::create_lcao_socket_solver(bool gamma_only)
{
    if (gamma_only)
    {
        return new SocketLCAO<double>();
    }
    return new SocketLCAO<std::complex<double>>();
}
#endif
