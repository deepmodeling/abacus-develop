#ifndef SYMMETRY_RHO_H
#define SYMMETRY_RHO_H
#include "module_elecstate/module_charge/charge.h"
#include "module_pw/pw_basis.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

#include "module_cell/module_symmetry/symmetry.h"

class Symmetry_rho
{
	public:
	Symmetry_rho();
	~Symmetry_rho();

	void begin(const int &spin_now, const Charge &CHR, const ModulePW::PW_Basis *pw, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const;

	private:

	void psymm(double *rho_part, const ModulePW::PW_Basis *pw, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const;

};

#endif
