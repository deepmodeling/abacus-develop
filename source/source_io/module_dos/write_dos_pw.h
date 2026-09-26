#ifndef WRITE_DOS_PW_H
#define WRITE_DOS_PW_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/fp_energy.h"

namespace ModuleIO
{
	/// @brief calculate density of states(DOS) for PW base
	void write_dos_pw(
			const UnitCell& ucell,
			const ModuleBase::matrix &ekb,
			const ModuleBase::matrix &wg,
			const K_Vectors& kv,
			const int nbands,
			const int istep_in,
			const elecstate::Efermi &energy_fermi,
			const double &dos_edelta_ev,
			const double &dos_scale,
			const double &bcoeff,
			const int nspin,
			const int out_dos,
			const bool dos_setemax,
			const double dos_emax_ev,
			const bool dos_setemin,
			const double dos_emin_ev,
			const bool two_fermi,
			const bool out_app_flag,
			const int bndpar,
			const std::string& global_out_dir,
			std::ofstream& ofs_running);
}
#endif
