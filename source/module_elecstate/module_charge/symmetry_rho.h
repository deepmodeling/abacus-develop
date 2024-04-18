#ifndef SYMMETRY_RHO_H
#define SYMMETRY_RHO_H
#include "module_elecstate/module_charge/charge.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_pw/hamilt_pwdft/parallel_grid.h"

#include "module_cell/module_symmetry/symmetry.h"

class Symmetry_rho
{
	public:
	Symmetry_rho();
	~Symmetry_rho();

    void begin(const int& spin_now,
               const Charge& CHR,
               const ModulePW::PW_Basis* rho_basis,
               Parallel_Grid& Pgrid,
               ModuleSymmetry::Symmetry& symm) const;

  private:
	// in real space:
	static void psymm(double *rho_part, const ModulePW::PW_Basis *rho_basis, Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) ;
	// in reciprocal space:
	void psymmg(std::complex<double>* rhog_part, const ModulePW::PW_Basis *rho_basis, 
			Parallel_Grid &Pgrid, ModuleSymmetry::Symmetry &symm) const;
#ifdef __MPI
	static void reduce_to_fullrhog(const ModulePW::PW_Basis *rho_basis, std::complex<double>* rhogtot, 
			std::complex<double>* rhogin, int* ig2isztot, const int* ig2iszin, int max_npw) ;
	static void rhog_piece_to_all(const ModulePW::PW_Basis *rho_basis,
			std::complex<double>* rhogtot, std::complex<double>* rhog_part) ;
#endif
	static void get_ixyz2ipw(const ModulePW::PW_Basis *rho_basis, const int* ig2isztot, 
			const int* fftixy2is, int* ixyz2ipw) ;	//(ix, iy, iz) -> (ip, ig)
};

#endif
