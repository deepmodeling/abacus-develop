#ifndef FORCES_H
#define FORCES_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"

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

    void init(ModuleBase::matrix& force);

private:

    int nat;
	static double output_acc;

    void cal_force_loc(ModuleBase::matrix& forcelc);
    void cal_force_ew(ModuleBase::matrix& forceion);
    void cal_force_cc(ModuleBase::matrix& forcecc);
    void cal_force_nl(ModuleBase::matrix& forcenl);
    void cal_force_scc(ModuleBase::matrix& forcescc);

    static void print( const std::string &name, const ModuleBase::matrix &f, bool rv=true );
    static void print_to_files( std::ofstream &ofs, const std::string &name, const ModuleBase::matrix &f );
};

#endif
