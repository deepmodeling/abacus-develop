#include "exx_abfs-parallel-communicate-dm3.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef __MPI
const ModuleBase::matrix &Exx_Abfs::Parallel::Communicate::DM3::D_phase(
	const ModuleBase::matrix &DK, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	assert(box2 == Abfs::Vector3_Order<int>(0,0,0));
	return DK;
}
ModuleBase::matrix Exx_Abfs::Parallel::Communicate::DM3::D_phase(
	const ModuleBase::ComplexMatrix &DK, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	return (DK * exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (box2*GlobalC::ucell.latvec)) )).real();
}

template<typename Tmatrix> std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
#endif
