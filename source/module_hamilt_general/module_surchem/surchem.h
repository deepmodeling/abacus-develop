#ifndef SURCHEM_H
#define SURCHEM_H

#include "../../module_base/atom_in.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"

class surchem
{
  public:
    surchem();
    ~surchem();

    double* TOTN_real;
    double* delta_phi;
    double* epspot;
    ModuleBase::matrix Vcav;
    ModuleBase::matrix Vel;
    double qs;

    // get atom info
    atom_in GetAtom;

    // allocate memory and deallocate them
    void allocate(const int& nrxx, const int& nspin);

    void clear();

    void cal_epsilon(const ModulePW::PW_Basis* rho_basis, const double* PS_TOTN_real, double* epsilon, double* epsilon0);

    void cal_pseudo(const UnitCell& cell,
                    const ModulePW::PW_Basis* rho_basis,
                    const complex<double>* Porter_g,
                    complex<double>* PS_TOTN,
                    Structure_Factor* sf);

    void gauss_charge(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, complex<double>* N, Structure_Factor* sf);

    void cal_totn(const UnitCell& cell,
                  const ModulePW::PW_Basis* rho_basis,
                  const complex<double>* Porter_g,
                  complex<double>* N,
                  complex<double>* TOTN,
                  const double* vlocal);

    void createcavity(const UnitCell& ucell,
                      const ModulePW::PW_Basis* rho_basis,
                      const complex<double>* PS_TOTN,
                      double* vwork);

    ModuleBase::matrix cal_vcav(const UnitCell& ucell,
                                const ModulePW::PW_Basis* rho_basis,
                                complex<double>* PS_TOTN,
                                int nspin);

    ModuleBase::matrix cal_vel(const UnitCell& cell,
                               const ModulePW::PW_Basis* rho_basis,
                               complex<double>* TOTN,
                               complex<double>* PS_TOTN,
                               int nspin);

    double cal_Ael(const UnitCell& cell,
                   const int& nrxx,  // num. of real space grids on current core
                   const int& nxyz); // total num. of real space grids

    double cal_Acav(const UnitCell& cell,
                    const int& nxyz); // total num. of real space grids

    void cal_Acomp(const UnitCell& cell,
                   const ModulePW::PW_Basis* rho_basis,
                   const double* const* const rho,
                   std::vector<double>& res);

    void minimize_cg(const UnitCell& ucell,
                     const ModulePW::PW_Basis* rho_basis,
                     double* d_eps,
                     const complex<double>* tot_N,
                     complex<double>* phi,
                     int& ncgsol);

    void Leps2(const UnitCell& ucell,
               const ModulePW::PW_Basis* rho_basis,
               complex<double>* phi,
               double* epsilon,            // epsilon from shapefunc, dim=nrxx
               complex<double>* gradphi_x, // dim=ngmc
               complex<double>* gradphi_y,
               complex<double>* gradphi_z,
               complex<double>* phi_work,
               complex<double>* lp);

    ModuleBase::matrix v_correction(const UnitCell& cell,
                                    const ModulePW::PW_Basis* rho_basis,
                                    const int& nspin,
                                    const double* const* const rho,
                                    const double* vlocal,
                                    Structure_Factor* sf);

    void test_V_to_N(ModuleBase::matrix& v,
                     const UnitCell& cell,
                     const ModulePW::PW_Basis* rho_basis,
                     const double* const* const rho);

    void cal_force_sol(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, ModuleBase::matrix& forcesol);

    void get_totn_reci(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, complex<double>* totn_reci);

    void induced_charge(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, double* induced_rho);

  private:
};

namespace GlobalC
{
extern surchem solvent_model;
}

#endif
