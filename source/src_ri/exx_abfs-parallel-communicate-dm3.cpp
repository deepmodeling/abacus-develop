#include "exx_abfs-parallel-communicate-dm3.h"
#include "exx_abfs-parallel-communicate-dm3-template.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"
#include "abfs-template.h"

#ifdef __MPI
void Exx_Abfs::Parallel::Communicate::DM3::cal_DM(const double threshold_D,
    Local_Orbital_Charge &loc)
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM3::cal_DM");
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DR_a2D = GlobalV::GAMMA_ONLY_LOCAL
		? K_to_R(loc.dm_gamma, threshold_D, *loc.ParaV)
		: K_to_R(loc.dm_k, threshold_D, *loc.ParaV);
	DMr = allreduce.a2D_to_exx(DR_a2D);

}
#endif
