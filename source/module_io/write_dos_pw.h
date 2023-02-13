#ifndef WRITE_DOS_PW_H
#define WRITE_DOS_PW_H

#include "module_elecstate/elecstate.h"
namespace ModuleIO
{
	/// @brief calculate density of states(DOS) for PW base
	void write_dos_pw(const elecstate::ElecState* pelec,
		const double &dos_edelta_ev,
		const double &dos_scale,
		const double &ef);
}
#endif
