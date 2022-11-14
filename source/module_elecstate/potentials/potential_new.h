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
    //interface for SCF-converged, etxc vtxc for Energy, vnew for force_scc 
    void get_vnew(const Charge* chg, ModuleBase::matrix& vnew);
    
    //interfaces to get values
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
    double* get_fixed_v()
    {
        return this->v_effective_fixed.data();
    }
    const double* get_fixed_v() const
    {
        return this->v_effective_fixed.data();
    }

    //interface for printing
    // mohan add 2011-02-28
    // here vh is std::complex because the array is got after std::complex FFT.
    void write_potential(const int &is,
                         const int &iter,
                         const std::string &fn,
                         const ModuleBase::matrix &v,
                         const int &precision,
                         const int &hartree = 0) const;

    void write_elecstat_pot(const std::string &fn, const std::string &fn_ave, ModulePW::PW_Basis *rho_basis);

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