#ifndef SURCHEM_H
#define SURCHEM_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_cell/unitcell.h"
#include "../module_pw/pw_basis.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_pw/global.h"
#include "../src_pw/structure_factor.h"
#include "../src_pw/use_fft.h"
#include "atom_in.h"

class surchem
{
  public:
    surchem();
    ~surchem();

    double *TOTN_real;
    double *delta_phi;
    double *epspot;
    ModuleBase::matrix Vcav;
    ModuleBase::matrix Vel;
    double qs;

    // compensating charge (in real space, used to cal_Acomp)
    double *comp_real;
    double *phi_comp_R;

    // compensating charge params
    double comp_q;
    double comp_l;
    double comp_center;
    int comp_dim;
    // get atom info
    atom_in GetAtom;

    void allocate(const int &nrxx, const int &nspin);

    void cal_epsilon(ModulePW::PW_Basis *rho_basis, const double *PS_TOTN_real, double *epsilon, double *epsilon0);

    void cal_pseudo(const UnitCell &cell,
                    ModulePW::PW_Basis *rho_basis,
                    const complex<double> *Porter_g,
                    complex<double> *PS_TOTN);

    void add_comp_chg(const UnitCell &cell,
                      ModulePW::PW_Basis *rho_basis,
                      double q,
                      double l,
                      double center,
                      complex<double> *NG,
                      int dim,
                      bool flag); // Set value of comp_reci[ig_gge0] when flag is true. 

    void cal_comp_force(ModuleBase::matrix &force_comp, ModulePW::PW_Basis *rho_basis);

    void gauss_charge(const UnitCell &cell, ModulePW::PW_Basis *rho_basis, complex<double> *N);

    void cal_totn(const UnitCell &cell,
                  ModulePW::PW_Basis *rho_basis,
                  const complex<double> *Porter_g,
                  complex<double> *N,
                  complex<double> *TOTN);
    void createcavity(const UnitCell &ucell,
                      ModulePW::PW_Basis *rho_basis,
                      const complex<double> *PS_TOTN,
                      double *vwork);

    ModuleBase::matrix cal_vcav(const UnitCell &ucell,
                                ModulePW::PW_Basis *rho_basis,
                                complex<double> *PS_TOTN,
                                int nspin);

    ModuleBase::matrix cal_vel(const UnitCell &cell,
                               ModulePW::PW_Basis *rho_basis,
                               complex<double> *TOTN,
                               complex<double> *PS_TOTN,
                               int nspin);

    double cal_Ael(const UnitCell &cell, ModulePW::PW_Basis *rho_basis);

    double cal_Acav(const UnitCell &cell, ModulePW::PW_Basis *rho_basis);

    void cal_Acomp(const UnitCell &cell,
                   ModulePW::PW_Basis *rho_basis,
                   const double *const *const rho,
                   vector<double> &res);

    void minimize_cg(const UnitCell &ucell,
                     ModulePW::PW_Basis *rho_basis,
                     double *d_eps,
                     const complex<double> *tot_N,
                     complex<double> *phi,
                     int &ncgsol);

    void Leps2(const UnitCell &ucell,
               ModulePW::PW_Basis *rho_basis,
               complex<double> *phi,
               double *epsilon, // epsilon from shapefunc, dim=nrxx
               complex<double> *gradphi_x, // dim=ngmc
               complex<double> *gradphi_y,
               complex<double> *gradphi_z,
               complex<double> *phi_work,
               complex<double> *lp);

    ModuleBase::matrix v_correction(const UnitCell &cell,
                                    ModulePW::PW_Basis *rho_basis,
                                    const int &nspin,
                                    const double *const *const rho);

    ModuleBase::matrix v_compensating(const UnitCell &cell,
                                      ModulePW::PW_Basis *rho_basis,
                                      const int &nspin,
                                      const double *const *const rho);

    void test_V_to_N(ModuleBase::matrix &v,
                     const UnitCell &cell,
                     ModulePW::PW_Basis *rho_basis,
                     const double *const *const rho);

    void cal_force_sol(const UnitCell &cell, ModulePW::PW_Basis *rho_basis, ModuleBase::matrix &forcesol);

    void get_totn_reci(const UnitCell &cell, ModulePW::PW_Basis *rho_basis, complex<double> *totn_reci);

  private:
};

namespace GlobalC
{
extern surchem solvent_model;
}

#endif
