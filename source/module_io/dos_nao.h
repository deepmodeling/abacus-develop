#ifndef OUT_DOS_NAO_H
#define OUT_DOS_NAO_H

#include "module_io/nscf_fermi_surf.h"
#include "module_io/write_dos_lcao.h"
#include "module_elecstate/fp_energy.h"
#include "module_hamilt_general/hamilt.h"

namespace ModuleIO
{
	template<typename T>
	void out_dos_nao(
			const psi::Psi<T>* psi,
			hamilt::Hamilt<T>* p_ham,
			const Parallel_Orbitals &pv,
			const UnitCell& ucell,
			const K_Vectors& kv,
			const int nbands,
			const elecstate::efermi& eferm,
			const ModuleBase::matrix& ekb,
			const ModuleBase::matrix& wg,
			const double& dos_edelta_ev,
			const double& dos_scale,
			const double& dos_sigma);
}

#endif
