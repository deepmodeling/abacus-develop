#ifndef HHARTREEPW_H
#define HHARTREEPW_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"

namespace elecstate
{

class H_Hartree_pw
{
  public:
    H_Hartree_pw();
    ~H_Hartree_pw();

    // the Hartree energy
    static double hartree_energy;

    // compute the Hartree energy
    static ModuleBase::matrix v_hartree(const UnitCell &cell,
                                        ModulePW::PW_Basis *rho_basis,
                                        const int &nspin,
                                        const double *const *const rho);

    static int get_Z(string str);

    static void cast_C2R(complex<double> *src, double *dst, int dim);

    static void lapl_rho(const std::complex<double> *rhog, double *lapn);

    static void shape_gradn(const complex<double> *PS_TOTN, ModulePW::PW_Basis *rho_basis, double *eprime);

    static void eps_pot(const complex<double>* PS_TOTN,
                        const complex<double>* phi,
                        const ModulePW::PW_Basis* rho_basis,
                        double* d_eps,
                        double* vwork);

    static void test_res(const UnitCell &ucell,
                         ModulePW::PW_Basis *rho_basis,
                         const complex<double> *tot_N,
                         complex<double> *phi,
                         double *d_eps);

  private:
};

} // namespace elecstate

#include "pot_base.h"
namespace elecstate
{
// new interface for elecstate::Potential
class PotHartree : public PotBase
{
  public:
    PotHartree(const ModulePW::PW_Basis* rho_basis_in);

    void cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff);
};

} // namespace elecstate

#endif // Hartree energy
