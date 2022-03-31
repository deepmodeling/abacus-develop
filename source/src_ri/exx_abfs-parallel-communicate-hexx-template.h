#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H

#include "exx_abfs-parallel-communicate-hexx.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"
#include <cassert>
#ifdef __MPI
/*
template< typename T >
T Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D( const std::map<size_t,std::map<size_t,T>> & H_a2D ) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D");
	
	T H_m2D;s
	if(GlobalV::KS_SOLVER=="genelpa")
		H_m2D.create( GlobalC::ParaO.ncol, GlobalC::ParaO.nrow );
	else
		H_m2D.create( GlobalC::ParaO.nrow, GlobalC::ParaO.ncol );
	
	for( const auto H_a2DA : H_a2D )
	{
		const size_t iat1 = H_a2DA.first;
		for( const auto H_a2DB : H_a2DA.second )
		{
			const size_t iat2 = H_a2DB.first;
			const T & H = H_a2DB.second;
			
			for( int iw1=0; iw1!=H.nr; ++iw1 )
			{
				const int iwt1 = GlobalC::ucell.itiaiw2iwt( GlobalC::ucell.iat2it[iat1], GlobalC::ucell.iat2ia[iat1], iw1 );
				const int iwt1_m2D = GlobalC::ParaO.trace_loc_row[iwt1];
				if( iwt1_m2D == -1 )	continue;
				
				for( int iw2=0; iw2!=H.nc; ++iw2 )
				{
					const int iwt2 = GlobalC::ucell.itiaiw2iwt( GlobalC::ucell.iat2it[iat2], GlobalC::ucell.iat2ia[iat2], iw2 );
					const int iwt2_m2D = GlobalC::ParaO.trace_loc_col[iwt2];
					if( iwt2_m2D == -1 )	continue;
					
					if(GlobalV::KS_SOLVER=="genelpa")
						H_m2D( iwt2_m2D, iwt1_m2D ) = H(iw1,iw2);
					else
						H_m2D( iwt1_m2D, iwt2_m2D ) = H(iw1,iw2);
				}
			}
		}
	}
	return H_m2D;
}
*/

template<>
ModuleBase::matrix Exx_Abfs::Parallel::Communicate::Hexx::H_phase<ModuleBase::matrix>(
	ModuleBase::matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	assert(box2 == Abfs::Vector3_Order<int>(0,0,0));
	return std::forward<ModuleBase::matrix>(HR);
}

template<>
ModuleBase::ComplexMatrix Exx_Abfs::Parallel::Communicate::Hexx::H_phase<ModuleBase::ComplexMatrix>(
	ModuleBase::matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	return ModuleBase::ComplexMatrix(HR) * exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (box2*GlobalC::ucell.latvec)) );
}

template<typename Tmatrix>
Tmatrix Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D(const Parallel_Orbitals &pv, 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_a2D, const int ik) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D");

	Tmatrix HK_m2D;
	if(GlobalV::KS_SOLVER=="genelpa")
		HK_m2D.create( pv.ncol, pv.nrow );
	else
		HK_m2D.create( pv.nrow, pv.ncol );

	const int is_begin = (GlobalV::NSPIN==4) ? 0 : GlobalC::kv.isk[ik];
	const int is_end = (GlobalV::NSPIN==4) ? 4 : GlobalC::kv.isk[ik]+1;
	for(int is=is_begin; is!=is_end; ++is)
	{
		for(auto & HR_a2D_A : HR_a2D[is])
		{
			const size_t iat1 = HR_a2D_A.first;
			const size_t it1 = GlobalC::ucell.iat2it[iat1];
			const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
			for(auto & HR_a2D_B : HR_a2D_A.second)
			{
				const size_t iat2 = HR_a2D_B.first;
				const size_t it2 = GlobalC::ucell.iat2it[iat2];
				const size_t ia2 = GlobalC::ucell.iat2ia[iat2];

				Tmatrix HK_a2D;
				for(auto & HR_a2D_C : HR_a2D_B.second)
				{
					const Abfs::Vector3_Order<int> &box2 = HR_a2D_C.first;
					if(HK_a2D.c)
						HK_a2D += H_phase<Tmatrix>(std::forward<ModuleBase::matrix>(HR_a2D_C.second), ik, box2);
					else
						HK_a2D = H_phase<Tmatrix>(std::forward<ModuleBase::matrix>(HR_a2D_C.second), ik, box2);
				}

				for(int iw1=0; iw1!=HK_a2D.nr; ++iw1)
				{
					const int iw1_tmp = (GlobalV::NSPIN==4)
						? (iw1*2+is/2)
						: iw1;
					const int iwt1 = GlobalC::ucell.itiaiw2iwt(it1,ia1,iw1_tmp);
					const int iwt1_m2D = pv.trace_loc_row[iwt1];
					if(iwt1_m2D<0) continue;
					for(int iw2=0; iw2!=HK_a2D.nc; ++iw2)
					{
						const int iw2_tmp = (GlobalV::NSPIN==4)
							? (iw2*2+is%2)
							: iw2;
						const int iwt2 = GlobalC::ucell.itiaiw2iwt(it2,ia2,iw2_tmp);
						const int iwt2_m2D = pv.trace_loc_col[iwt2];
						if(iwt2_m2D<0) continue;

						if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
							HK_m2D(iwt2_m2D,iwt1_m2D) = HK_a2D(iw1,iw2);
						else
							HK_m2D(iwt1_m2D,iwt2_m2D) = HK_a2D(iw1,iw2);
					}
				}
			}
		}
	}
	return HK_m2D;
}

template<typename Tmatrix>
void Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D_mixing(const Parallel_Orbitals &pv, 
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> &HR_a2D,
	std::vector<Tmatrix> &HK_m2D,
	std::vector<std::deque<Tmatrix>> &HK_m2D_pulay_seq) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D_mixing");

	HK_m2D.resize(GlobalC::kv.nks);
	HK_m2D_pulay_seq.resize(GlobalC::kv.nks);
	for( int ik=0; ik!=GlobalC::kv.nks; ++ik )
	{
//gettimeofday( &t_start, NULL);
//			const std::map<size_t,std::map<size_t,ModuleBase::matrix>> HK_a2D = R_to_K(HR_a2D[ik]);
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::R_to_K\t"<<time_during(t_start)<<std::endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HK_a2D_"+ModuleBase::GlobalFunc::TO_STRING(ik)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), HK_a2D );
//gettimeofday( &t_start, NULL);
		switch(mixing_mode)
		{
			case Mixing_Mode::No:
				HK_m2D[ik] = Ra2D_to_Km2D<Tmatrix>(pv, HR_a2D, ik);
				break;
			case Mixing_Mode::Plain:
				if( HK_m2D[ik].nr && HK_m2D[ik].nc )
					HK_m2D[ik] = (1-mixing_beta) * HK_m2D[ik] + mixing_beta * Ra2D_to_Km2D<Tmatrix>(pv, HR_a2D, ik);
				else
					HK_m2D[ik] = Ra2D_to_Km2D<Tmatrix>(pv, HR_a2D, ik);
				break;
			case Mixing_Mode::Pulay:
				HK_m2D[ik] = pulay_mixing( HK_m2D[ik], HK_m2D_pulay_seq[ik], Ra2D_to_Km2D<Tmatrix>(pv, HR_a2D, ik) );
				break;
			default:
				throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;
		}
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D\t"<<time_during(t_start)<<std::endl;
	}
//ofs_matrixes( exx_lcao.test_dir+"test-HK_m2D_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK), HK_m2D );
}

template<typename Tmatrix>
Tmatrix Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing(
	const Tmatrix &H_pulay_old, std::deque<Tmatrix> &H_seq, const Tmatrix &H_new ) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing");

	Tmatrix H_pulay;
	if( 0==GlobalC::CHR.totstep )
	{
		H_seq.clear();
		H_pulay = H_new;
	}
	else
	{
		H_seq.push_back( 
			(1-GlobalC::CHR.mixing_beta) * H_pulay_old + GlobalC::CHR.mixing_beta * H_new );
		if( H_seq.size() > GlobalC::CHR.dstep+1 )
			H_seq.pop_front();
		
		if( 1==H_seq.size() )
		{
			H_pulay = H_seq[0];
		}
		else
		{
			auto pos_mod = [](const int i, const int N){ return (i%N+N)%N; };
			auto index = [&](const int i) -> int
			{
				const int alpha_size = H_seq.size()-1;
				const int alpha_begin = GlobalC::CHR.idstep - alpha_size;
				return pos_mod( alpha_begin+i, GlobalC::CHR.dstep );
			};
			H_pulay = (1+GlobalC::CHR.alpha[index(H_seq.size()-2)]) * H_seq[H_seq.size()-1];
			for( int i=1; i<=H_seq.size()-2; ++i )
				H_pulay += ( GlobalC::CHR.alpha[index(i-1)] - GlobalC::CHR.alpha[index(i)] ) * H_seq[i];
			H_pulay -= GlobalC::CHR.alpha[index(0)] * H_seq[0];
		}
	}
	return H_pulay;
}
#endif
#endif