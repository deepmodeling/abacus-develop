#ifndef SYMMETRY_RHO_H
#define SYMMETRY_RHO_H
#include "../src_pw/charge_broyden.h"
#include "../src_pw/pw_basis.h"
#include "../src_parallel/parallel_grid.h"

#include "symmetry.h"

class Symmetry_rho
{
	public:
	Symmetry_rho();
	~Symmetry_rho();
	
	void begin(const int &spin_now, const Charge_Broyden &CHR, const PW_Basis &pw, const Parallel_Grid &Pgrid, const Symmetry &symm) const;

	private:

	void psymm(double *rho_part, const PW_Basis &pw, const Parallel_Grid &Pgrid, const Symmetry &symm) const;
	
};

#endif
