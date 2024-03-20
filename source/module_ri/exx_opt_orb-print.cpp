#include "exx_opt_orb.h"
#include "../module_hamilt_pw/hamilt_pwdft/global.h"
#include "exx_abfs-jle.h"

void Exx_Opt_Orb::print_matrix(
    const K_Vectors &kv,
    const std::string& file_name,
	const std::vector<RI::Tensor<double>> &matrix_Q, 
	const std::vector<std::vector<RI::Tensor<double>>> &matrix_S,
	const RI::Tensor<double> &matrix_V,
	const size_t TA, const size_t IA, const size_t TB, const size_t IB,
	const ModuleBase::Element_Basis_Index::Range &range_jles, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_jles) const
{
	auto print_header = [&]( std::ofstream &ofs )
	{
		ofs << GlobalC::ucell.lat0 << std::endl;

		ofs << GlobalC::ucell.latvec.e11 << " " << GlobalC::ucell.latvec.e12 << " " << GlobalC::ucell.latvec.e13 << std::endl;
		ofs << GlobalC::ucell.latvec.e21 << " " << GlobalC::ucell.latvec.e22 << " " << GlobalC::ucell.latvec.e23 << std::endl;
		ofs << GlobalC::ucell.latvec.e31 << " " << GlobalC::ucell.latvec.e32 << " " << GlobalC::ucell.latvec.e33 << std::endl;
		
		if( TA==TB )
		{
			ofs << 1 << " ntype" << std::endl;
			ofs << GlobalC::ucell.atoms[TA].label << " label" << std::endl;
			if( IA==IB )
			{
				ofs << 1 << " na" << std::endl;
				ofs << GlobalC::ucell.atoms[TA].tau[IA].x << " " 
					<< GlobalC::ucell.atoms[TA].tau[IA].y << " " 
					<< GlobalC::ucell.atoms[TA].tau[IA].z << std::endl;
			}
			else
			{
				ofs << 2 << " na" << std::endl;
				ofs << GlobalC::ucell.atoms[TA].tau[IA].x << " " 
					<< GlobalC::ucell.atoms[TA].tau[IA].y << " "
					<< GlobalC::ucell.atoms[TA].tau[IA].z << std::endl;
				ofs << GlobalC::ucell.atoms[TB].tau[IB].x << " " 
					<< GlobalC::ucell.atoms[TB].tau[IB].y << " " 
					<< GlobalC::ucell.atoms[TB].tau[IB].z << std::endl;
			}
		}
		else
		{
			ofs << 2 << " ntype" << std::endl;
			ofs << GlobalC::ucell.atoms[TA].label << " label" << std::endl;
			ofs << 1 << " na" << std::endl;
			ofs << GlobalC::ucell.atoms[TA].tau[IA].x << " " 
				<< GlobalC::ucell.atoms[TA].tau[IA].y << " " 
				<< GlobalC::ucell.atoms[TA].tau[IA].z << std::endl;
			ofs << GlobalC::ucell.atoms[TB].label << " label" << std::endl;
			ofs << 1 << " na" << std::endl;
			ofs << GlobalC::ucell.atoms[TB].tau[IB].x << " " 
				<< GlobalC::ucell.atoms[TB].tau[IB].y << " " 
				<< GlobalC::ucell.atoms[TB].tau[IB].z << std::endl;
		}
		
		// ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc" << std::endl; // mohan add 2009-09-08

		// this parameter determine the total number of jlq.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc_jlq" << std::endl;//mohan modify 2009-09-08

		if(TA==TB)
			ofs << GlobalC::ORB.Phi[TA].getRcut() << " rcut_Jlq" << std::endl;
		else
			ofs << GlobalC::ORB.Phi[TA].getRcut() << " " << GlobalC::ORB.Phi[TB].getRcut() << " rcut_Jlq" << std::endl;

		// mohan add 'smooth' and 'smearing_sigma' 2009-08-28
		ofs << 0 << " smooth" << std::endl;
		ofs << 0 << " smearing_sigma" << std::endl;

		ofs << Exx_Abfs::Jle::tolerence << " tolerence" << std::endl;

		ofs << Exx_Abfs::Jle::Lmax << " lmax" << std::endl;

		ofs << kv.nkstot << " nks" << std::endl;
		assert( matrix_V.shape[0] == matrix_V.shape[1] );
		ofs	<< matrix_V.shape[0] << " nbands" << std::endl;
		
		auto cal_sum_M = [&range_jles](size_t T) -> size_t
		{
			size_t sum_M = 0;
			for( size_t L = 0; L!=range_jles[T].size(); ++L )
				sum_M += range_jles[T][L].M;
			return sum_M;
		};
		const size_t nwfc = (TA==TB && IA==IB) ? cal_sum_M(TA) : cal_sum_M(TA)+cal_sum_M(TB);
		ofs	<< nwfc << " nwfc" << std::endl;
		
		const size_t ecut_numberA = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * GlobalC::ORB.Phi[TA].getRcut() / ModuleBase::PI ); // Rydberg Unit
		const size_t ecut_numberB = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * GlobalC::ORB.Phi[TB].getRcut() / ModuleBase::PI ); // Rydberg Unit
		if(TA==TB)
			ofs	<< ecut_numberA << " ne" << std::endl;
		else
			ofs	<< ecut_numberA << " " << ecut_numberB << " ne" << std::endl;
		
		ofs << "<WEIGHT_OF_KPOINTS>" << std::endl;
		for( int ik=0; ik!=kv.nkstot; ++ik )		
		{
			ofs << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z;
			ofs << " " << kv.wk[ik] * 0.5 << std::endl;
		}
		ofs << "</WEIGHT_OF_KPOINTS>" << std::endl;

		ofs << std::endl;
	};
	
	
	auto print_Q = [&]( std::ofstream &ofs )
	{
		//---------------------
		//  < Psi | jY >
		//---------------------
		ofs<< "<OVERLAP_Q>" << std::endl;
		
		for( size_t ib=0; ib!=matrix_V.shape[0]; ++ib )
		{
			for( size_t iat=0; iat!=matrix_Q.size(); ++iat )
			{
				const size_t it = (iat==0) ? TA : TB;
				for( size_t il=0; il!=range_jles[it].size(); ++il )
				{
					for( size_t im=0; im!=range_jles[it][il].M; ++im )
					{
						for( size_t iq=0; iq!=range_jles[it][il].N; ++iq )
						{
							ofs<<matrix_Q[iat]( ib, index_jles[it][il][iq][im] )<<"\t"<<0<<std::endl;
						}
					}
				}
			}
		}
		
		ofs<< "</OVERLAP_Q>" << std::endl << std::endl;
	};
	
	
	auto print_S = [&]( std::ofstream &ofs, const double scale=1 )
	{
		//---------------------
		//  < jY | jY >
		//---------------------
		ofs<< "<OVERLAP_Sq>" <<std::endl;
		
		for( size_t iat1=0; iat1!=matrix_S.size(); ++iat1 )
		{
			const size_t it1 = (iat1==0) ? TA : TB;
			for( size_t il1=0; il1!=range_jles[it1].size(); ++il1 )
			{
				for( size_t im1=0; im1!=range_jles[it1][il1].M; ++im1 )
				{
					for( size_t iat2=0; iat2!=matrix_S[iat1].size(); ++iat2 )
					{
						const size_t it2 = (iat2==0) ? TA : TB;
						for( size_t il2=0; il2!=range_jles[it2].size(); ++il2 )
						{
							for( size_t im2=0; im2!=range_jles[it2][il2].M; ++im2 )
							{
								for( size_t iq1=0; iq1!=range_jles[it1][il1].N; ++iq1 )
								{
									for( size_t iq2=0; iq2!=range_jles[it2][il2].N; ++iq2 )
									{
										ofs<<matrix_S[iat1][iat2]( index_jles[it1][il1][iq1][im1], index_jles[it2][il2][iq2][im2] )*scale<<"\t"<<0<<std::endl;
									}
								}
							}
						}
					}
				}
			}
		}
		
		ofs<< "</OVERLAP_Sq>" << std::endl << std::endl;
	};
	
	
	auto print_V = [&]( std::ofstream &ofs, const double scale=1 )
	{
		//---------------------
		//  < Psi | Psi >
		//---------------------	
		ofs << "<OVERLAP_V>" << std::endl;
		
		for( size_t ib1=0; ib1!=matrix_V.shape[0]; ++ib1 )
		{
			for( size_t ib2=0; ib2!=matrix_V.shape[1]; ++ib2 )
			{
				ofs<<matrix_V(ib1,ib2)*scale<<"\t";
			}
			ofs<<std::endl;
		}
		
		ofs << "</OVERLAP_V>" << std::endl << std::endl;
	};
	
	std::ofstream ofs(file_name+"_"+ModuleBase::GlobalFunc::TO_STRING(TA)+"_"+ModuleBase::GlobalFunc::TO_STRING(IA)+"_"+ModuleBase::GlobalFunc::TO_STRING(TB)+"_"+ModuleBase::GlobalFunc::TO_STRING(IB));
	print_header(ofs);
	print_Q(ofs);
	print_S(ofs);
	print_V(ofs);
	ofs.close();
}