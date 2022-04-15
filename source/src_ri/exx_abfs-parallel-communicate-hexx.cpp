#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-parallel-communicate-hexx-template.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "exx_abfs-io.h"

#if EXX_H_COMM==1
#include "exx_abfs-parallel-communicate-allreduce-template.h"
#endif

#include "exx_abfs-io-template.h"

#include "../src_external/src_test/src_ri/exx_lcao-test.h"
#include "../src_external/src_test/test_function.h"
//#include <gperftools/profiler.h>

#ifdef __MPI
void Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D(const Parallel_Orbitals &pv, 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_exx,
	const std::pair<bool,bool> &io_HR_a2D )
{
	/*{
		static int istep=0;
		std::ofstream ofs("HR_exx_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<HR_exx<<std::endl;
	}*/

	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D");
	
//std::ofstream ofs_time("time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;	

//gettimeofday( &t_start, NULL);	
	MPI_Barrier(MPI_COMM_WORLD);		// Peize Lin test
//ofs_time<<"TIME@ MPI_Barrier\t"<<time_during(t_start)<<std::endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_exx_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), HR_exx );
	
//gettimeofday( &t_start, NULL);
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
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::exx_to_a2D\t"<<time_during(t_start)<<std::endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_a2D_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), HR_a2D );

	/*{
		static int istep=0;
		std::ofstream ofs("HR_a2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<HR_a2D<<std::endl;
	}*/

	if(GlobalV::GAMMA_ONLY_LOCAL)
		Ra2D_to_Km2D_mixing(pv, HR_a2D, HK_Gamma_m2D, HK_Gamma_m2D_pulay_seq);
	else
		Ra2D_to_Km2D_mixing(pv, HR_a2D, HK_K_m2D, HK_K_m2D_pulay_seq);

	/*{
		static int istep=0;
		std::ofstream ofs("HK_m2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		if(GlobalV::GAMMA_ONLY_LOCAL)
			ofs<<HK_Gamma_m2D<<std::endl;
		else
			ofs<<HK_K_m2D<<std::endl;
	}*/
//ofs_time.close();
}
#endif
/*
std::map<size_t,std::map<size_t,matrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,matrix>>> & HR) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	std::map<size_t,std::map<size_t,matrix>> HK;
	for( auto & HR1 : HR )
	{
		const size_t iat1 = HR1.first;
		for( auto & HR2 : HR1.second )
		{
			const size_t iat2 = HR2.first;
			assert(HR2.second.size()==1);
			HK[iat1][iat2] = std::move( HR2.second.at(ModuleBase::Vector3<int>{0,0,0}) );
		}
	}
	return HK;
}
*/

/*
std::map<size_t,std::map<size_t,ModuleBase::ComplexMatrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	const std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,matrix>>> & HR,
	const size_t ik) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	std::map<size_t,std::map<size_t,ModuleBase::ComplexMatrix>> HK;
	for( auto & HRA : HR )
	{
		const size_t iat1 = HRA.first;
		for( auto & HRB : HRA.second )
		{
			const size_t iat2 = HRB.first;
			
			ModuleBase::ComplexMatrix HK_tmp(
				GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw,
				GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
			for( auto & HRC : HRB.second )
			{
				const Abfs::Vector3_Order<int> & box2 = HRC.first;
				HK_tmp += ModuleBase::ComplexMatrix(HRC.second) * exp( ModuleBase::TWO_PI*IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (box2*GlobalC::ucell.latvec)) );
			}
			HK[iat1][iat2] = std::move(HK_tmp);
		}
	}
	return HK;
}
*/

