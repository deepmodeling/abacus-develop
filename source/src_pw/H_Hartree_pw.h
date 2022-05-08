#ifndef H_HARTREE_PW_H
#define H_HARTREE_PW_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_cell/unitcell.h"
#include "pw_basis.h"
#include "use_fft.h"

class H_Hartree_pw
{
  public:
    H_Hartree_pw();
    ~H_Hartree_pw();

    // the Hartree energy
    static double hartree_energy;

    // compute the Hartree energy
    static ModuleBase::matrix v_hartree(const UnitCell &cell,
                                        PW_Basis &pwb,
                                        const int &nspin,
                                        const double *const *const rho);

    static void gauss_charge(const UnitCell &cell, PW_Basis &pwb, complex<double> *N, const int flag);

    static int get_Z(string str);

    static ModuleBase::matrix v_correction(const UnitCell &cell,
                                           PW_Basis &pwb,
                                           const int &nspin,
                                           const double *const *const rho);

    static void cast_C2R(complex<double> *src, double *dst, int dim);

    static void Leps(const UnitCell &ucell,
                     PW_Basis &pwb,
                     complex<double> *phi,
                     double *epsilon, // epsilon from shapefunc
                     complex<double> *gradphi_x,
                     complex<double> *gradphi_y,
                     complex<double> *gradphi_z,
                     complex<double> *phi_work,
                     complex<double> *lp // output
    );

    static void Leps2(const UnitCell &ucell,
                      PW_Basis &pwb,
                      complex<double> *phi,
                      double *epsilon, // epsilon from shapefunc, dim=nrxx
                      complex<double> *gradphi_x, // dim=ngmc
                      complex<double> *gradphi_y,
                      complex<double> *gradphi_z,
                      complex<double> *phi_work,
                      complex<double> *lp);

    static void minimize(const UnitCell &ucell,
                         PW_Basis &pwb,
                         double *d_eps,
                         const complex<double> *tot_N,
                         complex<double> *phi,
                         int &ncgsol);

    // static void minimize_test(
    // 	const UnitCell &ucell,
    // 	PW_Basis &pwb,
    // 	double *d_eps,
    // 	const complex<double>* tot_N,
    // 	complex<double> *phi,
    // 	int &ncgsol);

    static void createcavity(const UnitCell &ucell,
                             PW_Basis &pwb,
                             const complex<double> *PS_TOTN,
                             double *vwork,
                             double &Acav);

    static void lapl_rho(const std::complex<double> *rhog, double *lapn);

    static void shape_gradn(const complex<double> *PS_TOTN, PW_Basis &pw, double *eprime);

    static void eps_pot(const complex<double> *PS_TOTN,
                        const complex<double> *phi,
                        PW_Basis &pw,
                        double *d_eps,
                        double *vwork);

    static void test_res(const UnitCell &ucell,
                         PW_Basis &pwb,
                         const complex<double> *tot_N,
                         complex<double> *phi,
                         double *d_eps);

  private:
};

#endif // Hartree energy
