#include "exx_opt_orb.h"
#include "../../source_pw/module_pwdft/global.h"
#include "exx_abfs-jle.h"
#include "source_base/tool_title.h"
#include <iomanip>

void Exx_Opt_Orb::print_matrix(
	const Exx_Info::Exx_Info_Opt_ABFs &info,
	const UnitCell& ucell,
	const K_Vectors &kv,
	const int Lmax,
	const std::vector<std::size_t> &ecut_number,
	const std::string &file_name,
	const std::vector<RI::Tensor<double>> &matrix_Q, 
	const std::vector<std::vector<RI::Tensor<double>>> &matrix_S,
	const RI::Tensor<double> &matrix_V,
	const std::size_t TA, const std::size_t IA, const std::size_t TB, const std::size_t IB,
	const std::vector<double>& orb_cutoff,
	const ModuleBase::Element_Basis_Index::Range &range_jles, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_jles) const
{
	auto print_header = [&]( std::ofstream &ofs )
	{
		ofs << ucell.lat0 << std::endl;

		ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << std::endl;
		ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << std::endl;
		ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << std::endl;
		
		if( TA==TB )
		{
			ofs << 1 << " ntype" << std::endl;
			ofs << ucell.atoms[TA].label << " label" << std::endl;
			if( IA==IB )
			{
				ofs << 1 << " na" << std::endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " " 
					<< ucell.atoms[TA].tau[IA].z << std::endl;
			}
			else
			{
				ofs << 2 << " na" << std::endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " "
					<< ucell.atoms[TA].tau[IA].z << std::endl;
				ofs << ucell.atoms[TB].tau[IB].x << " " 
					<< ucell.atoms[TB].tau[IB].y << " " 
					<< ucell.atoms[TB].tau[IB].z << std::endl;
			}
		}
		else
		{
			ofs << 2 << " ntype" << std::endl;
			ofs << ucell.atoms[TA].label << " label" << std::endl;
			ofs << 1 << " na" << std::endl;
			ofs << ucell.atoms[TA].tau[IA].x << " " 
				<< ucell.atoms[TA].tau[IA].y << " " 
				<< ucell.atoms[TA].tau[IA].z << std::endl;
			ofs << ucell.atoms[TB].label << " label" << std::endl;
			ofs << 1 << " na" << std::endl;
			ofs << ucell.atoms[TB].tau[IB].x << " " 
				<< ucell.atoms[TB].tau[IB].y << " " 
				<< ucell.atoms[TB].tau[IB].z << std::endl;
		}
		
		ofs << info.ecut_exx << " ecutwfc" << std::endl;

		// this parameter determine the total number of jlq.
		ofs << info.ecut_exx << " ecutwfc_jlq" << std::endl;

		if(TA==TB)
			{ ofs << orb_cutoff[TA] << " rcut_Jlq" << std::endl; }
		else
			{ ofs << orb_cutoff[TA] << " " << orb_cutoff[TB] << " rcut_Jlq" << std::endl; }

		ofs << 0 << " smooth" << std::endl;
		ofs << 0 << " smearing_sigma" << std::endl;

		ofs << info.tolerence << " tolerence" << std::endl;

		ofs << Lmax << " lmax" << std::endl;

		ofs << kv.get_nkstot() << " nks" << std::endl;
		assert( matrix_V.shape[0]*matrix_V.shape[1] == matrix_V.shape[2]*matrix_V.shape[3] );
		ofs	<< matrix_V.shape[0]*matrix_V.shape[1] << " nbands" << std::endl;
		
		auto cal_sum_M = [&range_jles](std::size_t T) -> std::size_t
		{
			std::size_t sum_M = 0;
			for( std::size_t L = 0; L!=range_jles[T].size(); ++L )
				{ sum_M += range_jles[T][L].M; }
			return sum_M;
		};
		const std::size_t nwfc = (TA==TB && IA==IB) ? cal_sum_M(TA) : cal_sum_M(TA)+cal_sum_M(TB);
		ofs	<< nwfc << " nwfc" << std::endl;
		
		for(const std::size_t ne : ecut_number)
			{ ofs << ne << " "; }
		ofs << "ne" << std::endl;
		
		ofs << "<WEIGHT_OF_KPOINTS>" << std::endl;
		for( int ik=0; ik!=kv.get_nkstot(); ++ik )		
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
		for( std::size_t iw0=0; iw0!=matrix_V.shape[0]; ++iw0 )
		{
			for( std::size_t iw1=0; iw1!=matrix_V.shape[1]; ++iw1 )
			{
				for( std::size_t iat=0; iat!=matrix_Q.size(); ++iat )
				{
					const std::size_t it = (iat==0) ? TA : TB;
					for( std::size_t il=0; il!=range_jles[it].size(); ++il )
					{
						for( std::size_t im=0; im!=range_jles[it][il].M; ++im )
						{
							for( std::size_t iq=0; iq!=range_jles[it][il].N; ++iq )
							{
								ofs<<matrix_Q[iat]( iw0, iw1, index_jles[it][il][iq][im] )<<"\t"<<0<<std::endl;
							}
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
		for( std::size_t iat1=0; iat1!=matrix_S.size(); ++iat1 )
		{
			const std::size_t it1 = (iat1==0) ? TA : TB;
			for( std::size_t il1=0; il1!=range_jles[it1].size(); ++il1 )
			{
				for( std::size_t im1=0; im1!=range_jles[it1][il1].M; ++im1 )
				{
					for( std::size_t iat2=0; iat2!=matrix_S[iat1].size(); ++iat2 )
					{
						const std::size_t it2 = (iat2==0) ? TA : TB;
						for( std::size_t il2=0; il2!=range_jles[it2].size(); ++il2 )
						{
							for( std::size_t im2=0; im2!=range_jles[it2][il2].M; ++im2 )
							{
								for( std::size_t iq1=0; iq1!=range_jles[it1][il1].N; ++iq1 )
								{
									for( std::size_t iq2=0; iq2!=range_jles[it2][il2].N; ++iq2 )
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
		for( std::size_t iw0=0; iw0!=matrix_V.shape[0]; ++iw0 )
		{
			for( std::size_t iw1=0; iw1!=matrix_V.shape[1]; ++iw1 )
			{
				for( std::size_t iw2=0; iw2!=matrix_V.shape[2]; ++iw2 )
				{
					for( std::size_t iw3=0; iw3!=matrix_V.shape[3]; ++iw3 )
					{
						ofs<<matrix_V(iw0,iw1,iw2,iw3)*scale<<"\t";
					}
				}
				ofs<<std::endl;
			}
		}
		ofs << "</OVERLAP_V>" << std::endl << std::endl;
	};
	
	ModuleBase::TITLE("Exx_Opt_Orb","print_matrix");
	std::ofstream ofs(file_name+"_"+std::to_string(TA)+"_"+std::to_string(IA)+"_"+std::to_string(TB)+"_"+std::to_string(IB));
	print_header(ofs);
	ofs<<std::setprecision(15);
	print_Q(ofs);
	print_S(ofs);
	print_V(ofs);
}
