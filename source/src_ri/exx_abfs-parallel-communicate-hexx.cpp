#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-parallel-communicate-hexx-template.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "exx_abfs-io.h"

#if EXX_H_COMM==1
#include "exx_abfs-parallel-communicate-allreduce-template.h"
#endif

#include "exx_abfs-io-template.h"

//#include <gperftools/profiler.h>

#ifdef __MPI
void Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D(const Parallel_Orbitals &pv, 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_exx,
	const std::pair<bool,bool> &io_HR_a2D )
{

	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D");

	MPI_Barrier(MPI_COMM_WORLD);		// Peize Lin test

	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> HR_a2D;
	if(io_HR_a2D.first)
		HR_a2D = Exx_Abfs::IO::input_binary<std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>>(
			GlobalV::global_out_dir+"HR_exx_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
	else
	{
		#if EXX_H_COMM==1
			Allreduce allreduce(MPI_COMM_WORLD,HR_exx);
			HR_a2D = allreduce.exx_to_a2D();
		#elif EXX_H_COMM==2
			HR_a2D = allreduce2.exx_to_a2D(HR_exx);
		#else
			#error "EXX_H_COMM"
		#endif
	}
	if(io_HR_a2D.second)
		Exx_Abfs::IO::output_binary( HR_a2D, GlobalV::global_out_dir+"HR_exx_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK) );

	if(GlobalV::GAMMA_ONLY_LOCAL)
		Ra2D_to_Km2D_mixing(pv, HR_a2D, HK_Gamma_m2D, HK_Gamma_m2D_pulay_seq);
	else
		Ra2D_to_Km2D_mixing(pv, HR_a2D, HK_K_m2D, HK_K_m2D_pulay_seq);
}
#endif

