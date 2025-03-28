#ifndef GATEFIELD_H
#define GATEFIELD_H

#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"

namespace elecstate
{
class Gatefield
{
  public:
    Gatefield();
    ~Gatefield();

    static void add_gatefield(double *vltot,
                              const UnitCell &cell,
                              const ModulePW::PW_Basis *rho_basis,
                              const bool &linear,
                              const bool &quadratic);

    static double mopopla(double &zgate, double z, bool flag);

    static void compute_force(const UnitCell &cell, ModuleBase::matrix &fgate);

    static double etotgatefield; // energy for gatefield
    static double rho_surface; // surface charge density of charged plate
    static double zgate; // position of charged plate
    static bool relax; //
    static bool block; // add a block potential or not
    static double block_down; // low bound of the block
    static double block_up; // high bound of the block
    static double block_height; // height of the block
};

} // namespace elecstate
#include "pot_base.h"
namespace elecstate
{
// new interface for elecstate::Potential
class PotGate : public PotBase
{
  public:
    PotGate(const ModulePW::PW_Basis *rho_basis_in, const UnitCell *ucell_in) : ucell_(ucell_in)
    {
        this->rho_basis_ = rho_basis_in;
        this->fixed_mode = true;
        this->dynamic_mode = false;
    };

    void cal_fixed_v(double *vl_pseudo) override
    {
        Gatefield::add_gatefield(vl_pseudo, *ucell_, this->rho_basis_, true, true);
    }

  private:
    const UnitCell *ucell_ = nullptr;
};

} // namespace elecstate

#endif