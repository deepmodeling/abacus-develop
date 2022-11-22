#ifndef FORCES_H
#define FORCES_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_pw/pw_basis.h"
#include "module_psi/psi.h"
#include "charge.h"

class Forces
{
public:
	friend class Force_Stress_LCAO;
    /* This routine is a driver routine which compute the forces
     * acting on the atoms, the complete forces in plane waves
     * is computed from 4 main parts
     * (1) cal_force_loc: contribution due to local potential.
     * (2) cal_foce_ew: contribution due to ewald potential.
     * (3) cal_force_cc: contributino due to NLCC.
     * (4) cal_nl: contribution due to the non-local pseudopotential.
     * (4) cal_scc: contributino due to incomplete SCF calculation.
     */
    Forces();
    ~Forces();

    void init(ModuleBase::matrix& force, const ModuleBase::matrix& wg, const Charge* const chr, const psi::Psi<std::complex<double>>* psi_in=nullptr);

protected:

    int nat;
	static double output_acc;

    void cal_force_loc(ModuleBase::matrix& forcelc, ModulePW::PW_Basis* rho_basis, const Charge* const chr);
    void cal_force_ew(ModuleBase::matrix& forceion, ModulePW::PW_Basis* rho_basis);
    void cal_force_cc(ModuleBase::matrix& forcecc, ModulePW::PW_Basis* rho_basis, const Charge* const chr);
    void cal_force_nl(ModuleBase::matrix& forcenl, const ModuleBase::matrix& wg, const psi::Psi<std::complex<double>>* psi_in=nullptr);
    void cal_force_scc(ModuleBase::matrix& forcescc, ModulePW::PW_Basis* rho_basis);

    static void print( const std::string &name, const ModuleBase::matrix &f, bool rv=true );
    static void print_to_files( std::ofstream &ofs, const std::string &name, const ModuleBase::matrix &f );
};

#endif
