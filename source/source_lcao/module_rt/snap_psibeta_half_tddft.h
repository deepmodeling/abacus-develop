#ifndef SNAP_PSIBETA_HALF_TDDFT
#define SNAP_PSIBETA_HALF_TDDFT

#include <vector>
#include <complex>

#include "source_base/vector3.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_cell/setup_nonlocal.h"

namespace module_rt
{
	// calculate the tddft nonlocal potential term
	void snap_psibeta_half_tddft(
		const LCAO_Orbitals &orb,
		const InfoNonlocal &infoNL_,
		std::vector<std::vector<std::complex<double>>> &nlm,
		const ModuleBase::Vector3<double> &R1,
		const int &T1,
		const int &L1,
		const int &m1,
		const int &N1,
		const ModuleBase::Vector3<double> &R0, // The projector.
		const int &T0,
		const ModuleBase::Vector3<double> &A,
		const bool &calc_r
    );

} // namespace module_rt

#endif
