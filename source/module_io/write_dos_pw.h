#ifndef WRITE_DOS_PW_H
#define WRITE_DOS_PW_H
#include "module_base/matrix.h"
#include "module_cell/klist.h"

namespace ModuleIO
{
	/// @brief calculate density of states(DOS) for PW base
	void write_dos_pw(const ModuleBase::matrix &ekb,
		const ModuleBase::matrix &wg,
		const K_Vectors& kv,
		const double &dos_edelta_ev,
		const double &dos_scale,
		const double &bcoeff);
}
#endif
