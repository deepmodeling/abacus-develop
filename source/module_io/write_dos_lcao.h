#ifndef WRITE_DOS_LCAO_H
#define WRITE_DOS_LCAO_H
#include "module_psi/psi.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"

namespace ModuleIO
{
	/// @brief calculate density of states(DOS) and partial density of states(PDOS) and mulliken charge for LCAO base
	void write_dos_lcao(const psi::Psi<double> *psid,
		const psi::Psi<std::complex<double>> *psi,
		LCAO_Hamilt &uhm,
		const elecstate::ElecState* pelec,
		const int &out_dos,
		const double &dos_edelta_ev,
		const double &bcoeff,
		const double &dos_scale,
		const double &ef,
		const double &ef_up,
		const double &ef_dw);
}
#endif
