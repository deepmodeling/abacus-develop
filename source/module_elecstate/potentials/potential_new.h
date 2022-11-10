#ifndef POTENTIALNEW_H
#define POTENTIALNEW_H

#include "pot_base.h"
#include <vector>
#include "src_pw/VNL_in_pw.h"
#include "module_base/complexmatrix.h"

namespace elecstate
{

class Potential : public PotBase
{
  public:
    //In constructor, size of every potential components should be allocated
    Potential(
        const ModulePW::PW_Basis* rho_basis_in,
        const UnitCell_pseudo* ucell_in,
        const ModuleBase::matrix* vloc_in,
        const ModuleBase::ComplexMatrix* structure_factors_in,
        double* etxc_in,
        double* vtxc_in
    );
    ~Potential();

    //initialize potential when SCF begin
    void init_pot(int istep, const Charge* chg);
    //initialize potential components before SCF
    void pot_register(std::vector<std::string>& components_list);
    //update potential from current charge
    void update_from_charge(const Charge* chg, const UnitCell_pseudo* ucell);
    //update potential for purpose of TDDFT
    void update_for_tddft(int istep);
    
    ModuleBase::matrix& get_effective_v()
    {
        return this->v_effective;
    }
    const ModuleBase::matrix& get_effective_v()const
    {
        return this->v_effective;
    }
    ModuleBase::matrix& get_effective_vofk()
    {
        return this->vofk_effective;
    }
    const ModuleBase::matrix& get_effective_vofk()const
    {
        return this->vofk_effective;
    }
    ModuleBase::matrix& get_vnew()
    {
        return this->vnew;
    }
    const ModuleBase::matrix& get_vnew()const
    {
        return this->vnew;
    }

  private:

    void cal_v_eff(
      const Charge* chg, 
      const UnitCell_pseudo* ucell, 
      ModuleBase::matrix& v_eff) override;
    void cal_fixed_v(double *vl_pseudo) override;

    void allocate();

    std::vector<double> v_effective_fixed;
    ModuleBase::matrix v_effective;

    ModuleBase::matrix vofk_effective;

    //vnew(nspin,ncxyz) : V_out - V_in, needed in Force for PW base
    ModuleBase::matrix vnew;

    bool fixed_done = false;

    //gather etxc and vtxc in Potential, will be used in ESolver 
    double* etxc_ = nullptr;
    double* vtxc_ = nullptr;

    std::vector<PotBase*> components;

    const UnitCell_pseudo* ucell_ = nullptr;
    const ModuleBase::matrix* vloc_ = nullptr;
    const ModuleBase::ComplexMatrix* structure_factors_ = nullptr;
};

}

#endif