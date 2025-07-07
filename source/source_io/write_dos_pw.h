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
			const elecstate::efermi &energy_fermi,
			const double &dos_edelta_ev,
			const double &dos_scale,
			const double &bcoeff,
			std::ofstream& ofs_running);
}
#endif
