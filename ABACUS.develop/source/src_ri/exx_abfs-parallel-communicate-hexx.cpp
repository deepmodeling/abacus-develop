#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-parallel-communicate-hexx-template.h"
#include "src_pw/global.h"
#include "../module_base/global_function.h"
#include "exx_abfs-io.h"

#if EXX_H_COMM==1
#include "exx_abfs-parallel-communicate-allreduce-template.h"
#endif

#include "exx_abfs-io-template.h"

#include "src_external/src_test/src_ri/exx_lcao-test.h"
#include "src_external/src_test/test_function.h"
//#include <gperftools/profiler.h>

void Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D( 
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_exx,
	const pair<bool,bool> &io_HR_a2D )
{
	/*{
		static int istep=0;
		ofstream ofs("HR_exx_"+TO_STRING(istep++)+"_"+TO_STRING(MY_RANK));
		ofs<<HR_exx<<endl;
	}*/

	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D");
	
//ofstream ofs_time("time_"+TO_STRING(MY_RANK),ofstream::app);
//timeval t_start;	

//gettimeofday( &t_start, NULL);	
	MPI_Barrier(MPI_COMM_WORLD);		// Peize Lin test
//ofs_time<<"TIME@ MPI_Barrier\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_exx_"+TO_STRING(MY_RANK), HR_exx );
	
//gettimeofday( &t_start, NULL);
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> HR_a2D;
	if(io_HR_a2D.first)
		HR_a2D = Exx_Abfs::IO::input_binary<vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>>(
			global_out_dir+"HR_exx_"+TO_STRING(MY_RANK));
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
		Exx_Abfs::IO::output_binary( HR_a2D, global_out_dir+"HR_exx_"+TO_STRING(MY_RANK) );
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::exx_to_a2D\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_a2D_"+TO_STRING(MY_RANK), HR_a2D );

	/*{
		static int istep=0;
		ofstream ofs("HR_a2D_"+TO_STRING(istep++)+"_"+TO_STRING(MY_RANK));
		ofs<<HR_a2D<<endl;
	}*/

	if(GAMMA_ONLY_LOCAL)
		Ra2D_to_Km2D_mixing(HR_a2D, HK_Gamma_m2D, HK_Gamma_m2D_pulay_seq);
	else
		Ra2D_to_Km2D_mixing(HR_a2D, HK_K_m2D, HK_K_m2D_pulay_seq);

	/*{
		static int istep=0;
		ofstream ofs("HK_m2D_"+TO_STRING(istep++)+"_"+TO_STRING(MY_RANK));
		if(GAMMA_ONLY_LOCAL)
			ofs<<HK_Gamma_m2D<<endl;
		else
			ofs<<HK_K_m2D<<endl;
	}*/
//ofs_time.close();
}

/*
map<size_t,map<size_t,matrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	map<size_t,map<size_t,matrix>> HK;
	for( auto & HR1 : HR )
	{
		const size_t iat1 = HR1.first;
		for( auto & HR2 : HR1.second )
		{
			const size_t iat2 = HR2.first;
			assert(HR2.second.size()==1);
			HK[iat1][iat2] = std::move( HR2.second.at(Vector3<int>{0,0,0}) );
		}
	}
	return HK;
}
*/

/*
map<size_t,map<size_t,ComplexMatrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR,
	const size_t ik) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	map<size_t,map<size_t,ComplexMatrix>> HK;
	for( auto & HRA : HR )
	{
		const size_t iat1 = HRA.first;
		for( auto & HRB : HRA.second )
		{
			const size_t iat2 = HRB.first;
			
			ComplexMatrix HK_tmp(
				ucell.atoms[ucell.iat2it[iat1]].nw,
				ucell.atoms[ucell.iat2it[iat2]].nw);
			for( auto & HRC : HRB.second )
			{
				const Abfs::Vector3_Order<int> & box2 = HRC.first;
				HK_tmp += ComplexMatrix(HRC.second) * exp( TWO_PI*IMAG_UNIT * (kv.kvec_c[ik] * (box2*ucell.latvec)) );
			}
			HK[iat1][iat2] = std::move(HK_tmp);
		}
	}
	return HK;
}
*/

