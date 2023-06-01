#include "potential_new.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"

#include "module_elecstate/elecstate_getters.h"

#include <map>

namespace elecstate
{
Potential::Potential(const ModulePW::PW_Basis* rho_basis_in,
                     const UnitCell* ucell_in,
                     const ModuleBase::matrix* vloc_in,
                     Structure_Factor* structure_factors_in,
                     double* etxc_in,
                     double* vtxc_in)
    : ucell_(ucell_in), vloc_(vloc_in), structure_factors_(structure_factors_in), etxc_(etxc_in), vtxc_(vtxc_in)
{
    this->rho_basis_ = rho_basis_in;
    this->fixed_mode = true;
    this->dynamic_mode = true;

    // allocate memory for Potential.
    this->allocate();
}

Potential::~Potential()
{
    if (this->components.size() > 0)
    {
        for (auto comp: this->components)
        {
            delete comp;
        }
        this->components.clear();
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            delmem_sd_op()(gpu_ctx, s_v_effective);
            delmem_sd_op()(gpu_ctx, s_vofk_effective);
        }
        else {
            delmem_dd_op()(gpu_ctx, d_v_effective);
            delmem_dd_op()(gpu_ctx, d_vofk_effective);
        }
    }
    else {
        if (GlobalV::precision_flag == "single") {
            delmem_sh_op()(cpu_ctx, s_v_effective);
            delmem_sh_op()(cpu_ctx, s_vofk_effective);
        }
    }
}

void Potential::pot_register(std::vector<std::string>& components_list)
{
    ModuleBase::TITLE("Potential", "pot_register");
    // delete old components first.
    if (this->components.size() > 0)
    {
        for (auto comp: this->components)
        {
            delete comp;
        }
        this->components.clear();
    }

    // register components
    //---------------------------
    // mapping for register
    //---------------------------
    for (auto comp: components_list)
    {
        PotBase* tmp = this->get_pot_type(comp);
        this->components.push_back(tmp);
        //        GlobalV::ofs_running << "Successful completion of Potential's registration : " << comp << std::endl;
    }

    // after register, reset fixed_done to false
    this->fixed_done = false;

    return;
}

void Potential::allocate()
{
    ModuleBase::TITLE("Potential", "allocate");
    int nrxx = this->rho_basis_->nrxx;
    if (nrxx == 0)
        return;

    this->v_effective_fixed.resize(nrxx);
    ModuleBase::Memory::record("Pot::veff_fix", sizeof(double) * nrxx);

    this->v_effective.create(GlobalV::NSPIN, nrxx);
    ModuleBase::Memory::record("Pot::veff", sizeof(double) * GlobalV::NSPIN * nrxx);

    if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
    {
        this->vofk_effective.create(GlobalV::NSPIN, nrxx);
        ModuleBase::Memory::record("Pot::vofk", sizeof(double) * GlobalV::NSPIN * nrxx);
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            resmem_sd_op()(gpu_ctx, s_v_effective, GlobalV::NSPIN * nrxx);
            resmem_sd_op()(gpu_ctx, s_vofk_effective, GlobalV::NSPIN * nrxx);
        }
        else {
            resmem_dd_op()(gpu_ctx, d_v_effective, GlobalV::NSPIN * nrxx);
            resmem_dd_op()(gpu_ctx, d_vofk_effective, GlobalV::NSPIN * nrxx);
        }
    }
    else {
        if (GlobalV::precision_flag == "single") {
            resmem_sh_op()(cpu_ctx, s_v_effective, GlobalV::NSPIN * nrxx, "POT::sveff");
            resmem_sh_op()(cpu_ctx, s_vofk_effective, GlobalV::NSPIN * nrxx, "POT::svofk");
        }
        else {
            this->d_v_effective = this->v_effective.c;
            this->d_vofk_effective = this->vofk_effective.c;
        }
        // There's no need to allocate memory for double precision pointers while in a CPU environment
    }
}

void Potential::update_from_charge(const Charge* chg, const UnitCell* ucell)
{
    ModuleBase::TITLE("Potential", "update_from_charge");
    ModuleBase::timer::tick("Potential", "update_from_charge");
    if (!this->fixed_done)
    {
        this->cal_fixed_v(this->v_effective_fixed.data());
        this->fixed_done = true;
    }

    this->cal_v_eff(chg, ucell, this->v_effective);

    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, s_v_effective, this->v_effective.c, this->v_effective.nr * this->v_effective.nc);
            castmem_d2s_h2d_op()(gpu_ctx, cpu_ctx, s_vofk_effective, this->vofk_effective.c, this->vofk_effective.nr * this->vofk_effective.nc);
        }
        else {
            syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, d_v_effective, this->v_effective.c, this->v_effective.nr * this->v_effective.nc);
            syncmem_d2d_h2d_op()(gpu_ctx, cpu_ctx, d_vofk_effective, this->vofk_effective.c, this->vofk_effective.nr * this->vofk_effective.nc);
        }
    }
    else {
        if (GlobalV::precision_flag == "single") {
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, s_v_effective, this->v_effective.c, this->v_effective.nr * this->v_effective.nc);
            castmem_d2s_h2h_op()(cpu_ctx, cpu_ctx, s_vofk_effective, this->vofk_effective.c, this->vofk_effective.nr * this->vofk_effective.nc);
        }
        // There's no need to synchronize memory for double precision pointers while in a CPU environment
    }

    ModuleBase::timer::tick("Potential", "update_from_charge");
}

void Potential::cal_fixed_v(double* vl_pseudo)
{
    ModuleBase::TITLE("Potential", "cal_fixed_v");
    ModuleBase::timer::tick("Potential", "cal_fixed_v");
    this->v_effective_fixed.assign(this->v_effective_fixed.size(), 0.0);
    for (size_t i = 0; i < this->components.size(); i++)
    {
        if (this->components[i]->fixed_mode)
        {
            this->components[i]->cal_fixed_v(vl_pseudo);
        }
    }

    ModuleBase::timer::tick("Potential", "cal_fixed_v");
}

void Potential::cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff)
{
    ModuleBase::TITLE("Potential", "cal_v_eff");
    int nspin_current = this->v_effective.nr;
    int nrxx = this->v_effective.nc;
    ModuleBase::timer::tick("Potential", "cal_v_eff");
    // first of all, set v_effective to zero.
    this->v_effective.zero_out();

    // add fixed potential components
    // nspin = 2, add fixed components for all
    // nspin = 4, add fixed components on first colomn
    for (int i = 0; i < nspin_current; i++)
    {
        if (i == 0 || nspin_current == 2)
        {
            ModuleBase::GlobalFunc::COPYARRAY(this->v_effective_fixed.data(), this->get_effective_v(i), nrxx);
        }
    }

    // cal effective by every components
    for (size_t i = 0; i < this->components.size(); i++)
    {
        if (this->components[i]->dynamic_mode)
        {
            this->components[i]->cal_v_eff(chg, ucell, v_eff);
        }
    }

    ModuleBase::timer::tick("Potential", "cal_v_eff");
}

void Potential::init_pot(int istep, const Charge* chg)
{
    ModuleBase::TITLE("Potential", "init_pot");
    ModuleBase::timer::tick("Potential", "init_pot");

    assert(istep >= 0);
    // fixed components only calculated in the beginning of SCF
    this->fixed_done = false;

    this->update_from_charge(chg, this->ucell_);

    // plots
    // figure::picture(this->vr_eff1,GlobalC::rhopw->nx,GlobalC::rhopw->ny,GlobalC::rhopw->nz);
    ModuleBase::timer::tick("Potential", "init_pot");
    return;
}

void Potential::get_vnew(const Charge* chg, ModuleBase::matrix& vnew)
{
    ModuleBase::TITLE("Potential", "get_vnew");
    vnew.create(this->v_effective.nr, this->v_effective.nc);
    vnew = this->v_effective;

    this->update_from_charge(chg, this->ucell_);
    //(used later for scf correction to the forces )
    for (int iter = 0; iter < vnew.nr * vnew.nc; ++iter)
    {
        vnew.c[iter] = this->v_effective.c[iter] - vnew.c[iter];
    }

    return;
}

template <>
float * Potential::get_v_effective_data()
{
    return this->v_effective.nc > 0 ? this->s_v_effective : nullptr;
}

template <>
double * Potential::get_v_effective_data()
{
    return this->v_effective.nc > 0 ? this->d_v_effective : nullptr;
}

template <>
float * Potential::get_vofk_effective_data()
{
    return this->vofk_effective.nc > 0 ? this->s_vofk_effective : nullptr;
}

template <>
double * Potential::get_vofk_effective_data()
{
    return this->vofk_effective.nc > 0 ? this->d_vofk_effective : nullptr;
}

} // namespace elecstate