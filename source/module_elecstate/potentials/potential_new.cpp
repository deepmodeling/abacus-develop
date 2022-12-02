#include "potential_new.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"

#include <map>
#ifdef __LCAO
#include "src_lcao/ELEC_evolve.h"
#endif

#include "H_Hartree_pw.h"
#include "efield.h"
#include "gatefield.h"
#include "pot_local.h"
#include "pot_surchem.hpp"
#include "pot_xc.h"

namespace elecstate
{
Potential::Potential(const ModulePW::PW_Basis* rho_basis_in,
                     const UnitCell* ucell_in,
                     const ModuleBase::matrix* vloc_in,
                     const ModuleBase::ComplexMatrix* structure_factors_in,
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
    #if (defined(__CUDA) || defined(__ROCM))
    if (GlobalV::device_flag == "gpu") {
        delmem_var_op()(this->gpu_ctx, d_v_effective);
    }
    #endif
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
    std::map<string, int> pot_register_map
        = {{"local", 1}, {"hartree", 2}, {"xc", 3}, {"surchem", 4}, {"efield", 5}, {"gatefield", 6}};
    for (auto comp: components_list)
    {
        PotBase* tmp = nullptr;
        int key = pot_register_map[comp];
        switch (key)
        {
        case 1: //"local"
            tmp = new PotLocal(this->vloc_, this->structure_factors_, this->rho_basis_);
            break;
        case 2: //"hartree"
            tmp = new PotHartree(this->rho_basis_);
            break;
        case 3: //"xc"
            tmp = new PotXC(this->rho_basis_, this->etxc_, this->vtxc_, &(this->vofk_effective));
            break;
        case 4: //"surchem"
            tmp = new PotSurChem(this->rho_basis_, this->v_effective_fixed.data(), &GlobalC::solvent_model);
            break;
        case 5: //"efield"
            tmp = new PotEfield(this->rho_basis_, this->ucell_, GlobalV::DIP_COR_FLAG);
            break;
        case 6: //"gatefield"
            tmp = new PotGate(this->rho_basis_, this->ucell_);
            break;
        default:
            ModuleBase::WARNING_QUIT("Potential::Init", "Please input correct component of potential!");
            break;
        }
        this->components.push_back(tmp);
        GlobalV::ofs_running << "Successful completion of Potential's registration : " << comp << std::endl;
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
    ModuleBase::Memory::record("Potential", "v_effective_fixed", nrxx, "double");

    this->v_effective.create(GlobalV::NSPIN, nrxx);
    ModuleBase::Memory::record("Potential", "vr_eff", GlobalV::NSPIN * nrxx, "double");

    if (XC_Functional::get_func_type() == 3 || XC_Functional::get_func_type() == 5)
    {
        this->vofk_effective.create(GlobalV::NSPIN, nrxx);
        ModuleBase::Memory::record("Potential", "vofk", GlobalV::NSPIN * nrxx, "double");
    }
    #if (defined(__CUDA) || defined(__ROCM))
    if (GlobalV::device_flag == "gpu") {
        resmem_var_op()(this->gpu_ctx, d_v_effective, GlobalV::NSPIN * nrxx);
    }
    #endif
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

    #if (defined(__CUDA) || defined(__ROCM))
    if (GlobalV::device_flag == "gpu") {
        syncmem_var_h2d_op()(this->gpu_ctx, this->cpu_ctx, d_v_effective, this->v_effective.c, this->v_effective.nr * this->v_effective.nc);
    }
    #endif

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

#ifdef __LCAO
    if (ELEC_evolve::td_vext != 0 && istep < ELEC_evolve::td_timescale)
    {
        this->update_for_tddft(istep);
    }
#endif

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

} // namespace elecstate