#ifndef EFIELD_H
#define EFIELD_H

#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_general/module_surchem/surchem.h"

namespace elecstate
{
class Efield
{
  public:
    Efield();
    ~Efield();

    static ModuleBase::matrix add_efield(const UnitCell& cell,
                                         const ModulePW::PW_Basis* rho_basis,
                                         const int& nspin,
                                         const double* const* const rho,
                                         surchem& solvent);

    static double cal_elec_dipole(const UnitCell& cell,
                                  const ModulePW::PW_Basis* rho_basis,
                                  const int& nspin,
                                  const double* const* const rho,
                                  const double& bmod);

    static double cal_ion_dipole(const UnitCell& cell, const double& bmod);

    static double cal_induced_dipole(const UnitCell& cell,
                                     const ModulePW::PW_Basis* rho_basis,
                                     surchem& solvent,
                                     const double& bmod);

    static double saw_function(const double &a, const double &b, const double &x);

    static void compute_force(const UnitCell &cell, ModuleBase::matrix &fdip);

    static void prepare(const UnitCell &cell, double &latvec, double &area);

    static void autoset(const UnitCell& cell);

    static double etotefield; // dipole energy
    static double tot_dipole; // total dipole
    static int efield_dir; // 0, 1, 2 denotes x, y, z direction for dipole correction
    static double efield_pos_max; // the maximum position of the saw function
    static double efield_pos_dec; // the decrease region length of the saw function
    static double efield_amp; // field amplitude (in a.u.) (1 a.u. = 51.44 10^10 V/m)
    static double bvec[3];
    static double bmod;
};

} // namespace elecstate

#include "pot_base.h"
namespace elecstate
{
// new interface for elecstate::Potential
class PotEfield : public PotBase
{
  public:
    PotEfield(const ModulePW::PW_Basis* rho_basis_in, const UnitCell* ucell_in, bool dipole) : ucell_(ucell_in)
    {
        this->rho_basis_ = rho_basis_in;
        if (!dipole)
        {
            this->fixed_mode = true;
            this->dynamic_mode = false;
        }
        else
        {
            this->fixed_mode = false;
            this->dynamic_mode = true;
        }
    };

    void cal_fixed_v(double *vl_pseudo) override
    {
        ModuleBase::matrix v_efield(GlobalV::NSPIN, rho_basis_->nrxx);
        v_efield = Efield::add_efield(*ucell_,
                                      const_cast<const ModulePW::PW_Basis*>(rho_basis_),
                                      GlobalV::NSPIN,
                                      nullptr,
                                      GlobalC::solvent_model);
        for (int ir = 0; ir < rho_basis_->nrxx; ++ir)
        {
            vl_pseudo[ir] += v_efield(0, ir);
        }
    }

    void cal_v_eff(const Charge *chg, const UnitCell *ucell, ModuleBase::matrix &v_eff) override
    {
        v_eff += Efield::add_efield(*ucell,
                                    const_cast<const ModulePW::PW_Basis*>(rho_basis_),
                                    v_eff.nr,
                                    chg->rho,
                                    GlobalC::solvent_model);
    }

  private:
    const UnitCell *ucell_ = nullptr;
};

} // namespace elecstate

#endif