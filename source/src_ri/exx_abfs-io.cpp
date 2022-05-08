#include <fstream>
#include <stdexcept>

#include "exx_abfs-io.h"
#include "exx_abfs-jle.h"
#include "exx_abfs-abfs_index.h"
#include "../src_pw/global.h"
#include "../module_orbital/ORB_read.h"
#include "../module_base/global_function.h"
#include "../module_base/math_integral.h" // mohan add 2021-04-03


std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs(
	const LCAO_Orbitals &orbs,
	const std::vector<std::string> &files_abfs,
	const double kmesh_times )
{
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs( files_abfs.size() );
	for( size_t T=0; T!=files_abfs.size(); ++T )
		abfs[T] = construct_abfs_T( 
			files_abfs[T],
			T,
			static_cast<int>(orbs.get_kmesh() * kmesh_times) | 1,			// Nk must be odd
//			orbs.get_dk() / kmesh_times,
			orbs.get_dk(),								// Peize Lin change 2017-04-16
			orbs.get_dr_uniform() );		
	
	return abfs;
}

std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> Exx_Abfs::IO::construct_abfs( 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> & abfs_pre,
	const LCAO_Orbitals &orbs,
	const std::vector<std::string> &files_abfs,
	const double kmesh_times )
{
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> 
		&&abfs = construct_abfs( orbs, files_abfs, kmesh_times );
			
	assert( abfs.size() == abfs_pre.size() );
	for( size_t T=0; T!=abfs.size(); ++T )
	{
		if( abfs[T].size() < abfs_pre[T].size() )
			abfs[T].resize( abfs_pre[T].size() );
		for( size_t L=0; L!=abfs_pre[T].size(); ++L )
		{
			abfs[T][L].insert( abfs[T][L].begin(), abfs_pre[T][L].begin(), abfs_pre[T][L].end() );
		}
	}	
	
	return abfs;
}

std::vector<std::vector<Numerical_Orbital_Lm>> Exx_Abfs::IO::construct_abfs_T( 
	const std::string & file_name,
	const int &T,
	const int &nk,
	const double &dk,
	const double &dr_uniform)
{
	std::string label;
	size_t L_size;
	std::map<size_t,size_t> N_size;
	size_t meshr;
	double dr;
	std::map<size_t,std::map<size_t,std::vector<double>>> psis;
	
	/*----------------------
	  1.read abfs
	----------------------*/
	std::string word;
	
	std::ifstream ifs( file_name.c_str() );
	if(!ifs)
		throw std::runtime_error(" Can't find the abfs ORBITAL file.");	
	
	while( ifs.good() )
	{
		ifs >> word;
		
		if( "Element"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, label );
		}
		else if ( "Lmax"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, L_size );
			
			if( L_size>=9 )
			{
				std::stringstream ss;
				ss<<"Lmax>=9 error in "<<__FILE__<<" line "<<__LINE__;
				throw std::invalid_argument(ss.str());
			}
		}
		else if ( "Sorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[0] );
		}
		else if ( "Porbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[1] );
		}
		else if ( "Dorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[2] );
		}
		else if ( "Forbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[3] );
		}
		else if ( "Gorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[4] );
		}
		else if ( "Horbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[5] );
		}
		else if ( "Iorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[6] );
		}
		else if ( "Jorbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[7] );
		}
		else if ( "Korbital-->"==word )
		{
			ModuleBase::GlobalFunc::READ_VALUE( ifs, N_size[8] );
		}
		else if ( "END"==word )
		{
			break;
		}		
	}
	
	ModuleBase::CHECK_NAME(ifs, "Mesh");
	ifs >> meshr;
	
	ModuleBase::CHECK_NAME(ifs, "dr");
	ifs >> dr;

	while(ifs.good())
	{
		ifs >> word;
		if(word=="Type")
		{
			std::string s_L, s_N;
			ifs >> s_L >> s_N;
			
			size_t T,L,N;
			ifs >> T >> L >> N;
			
			psis[L][N].resize(meshr);
			for(int ir=0; ir!=meshr; ir++)
			{
				ifs >> psis[L][N][ir];
			}		
		}

	}
	ifs.close();

	 
	/*----------------------
	  2.check L,N orbital
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
		if( N_size.find(L) == N_size.end() )
		{
			std::stringstream ss;
			ss<<"Can't find N of L="<<L<<" in "<<file_name;
			throw std::domain_error(ss.str());
		}
	for( size_t L=0; L<=L_size; ++L )
		for( size_t N=0; N!=N_size[L]; ++N )
			if( psis.find(L)==psis.end() || psis[L].find(N)==psis[L].end() )
			{
				std::stringstream ss;
				ss<<"Can't find abf of L="<<L<<" T="<<T<<" in "<<file_name;
				throw std::domain_error(ss.str());
			}

			
	/*----------------------
	  3.rab, radial
	----------------------*/		
	if(meshr%2==0)	++meshr;
	
	std::vector<double> rab(meshr);
	std::vector<double> radial(meshr);
	for( int ir=0; ir!=meshr; ++ir )
	{
		rab[ir] = dr;
		radial[ir] = ir*dr; //mohan 2010-04-19
	}

	
	/*----------------------
	  4.normalize psi
	----------------------*/
	for( size_t L=0; L<=L_size; ++L )
	{	
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			std::vector<double> psir(meshr);
			std::vector<double> inner(meshr);
			psis[L][N].resize(meshr);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psir[ir] = psis[L][N][ir] * radial[ir];
				inner[ir] = psir[ir] * psir[ir];
			}
			double unit = 0.0;	
			ModuleBase::Integral::Simpson_Integral(meshr, ModuleBase::GlobalFunc::VECTOR_TO_PTR(inner), ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab), unit);
			for( int ir=0; ir!=meshr; ++ir )
			{
				psis[L][N][ir] /= sqrt(unit);
			}
		}
	}

	
	/*----------------------
	  5.construct abfs
	----------------------*/	
	std::vector<std::vector<Numerical_Orbital_Lm>> abfs_T;

	abfs_T.resize(L_size+1);
	for( size_t L=0; L<=L_size; ++L )
	{
		abfs_T[L].resize(N_size[L]);
		for( size_t N=0; N!=N_size[L]; ++N )
		{
			abfs_T[L][N].set_orbital_info(
				label,
				T, //type
				L, //angular momentum L
				N, // number of orbitals of this L
				meshr, // number of radial mesh
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(radial),// radial mesh value(a.u.)
				Numerical_Orbital_Lm::Psi_Type::Psi,
				ModuleBase::GlobalFunc::VECTOR_TO_PTR(psis[L][N]), // radial wave function
				nk,
				dk,
				dr_uniform,
				false,
				true, GlobalV::CAL_FORCE);		
		}
	}
	
	return abfs_T;
}

void Exx_Abfs::IO::print_matrix( 
		const std::string &file_name_prefix, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> &matrixes_Q, 
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrixes_S,
		const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> &matrixes_V,
		const ModuleBase::Element_Basis_Index::Range &range_jles, 
		const ModuleBase::Element_Basis_Index::IndexLNM &index_jles, 
		const ModuleBase::Element_Basis_Index::Range &range_lcaos,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos )
{
	auto print_header = [&]( std::ofstream &ofs, size_t TA, size_t IA, size_t TB, size_t IB )
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

		assert( GlobalC::ORB.Phi[TA].getRcut() == GlobalC::ORB.Phi[TB].getRcut() );
		ofs << GlobalC::ORB.Phi[TA].getRcut() << " rcut_Jlq" << std::endl;

		// mohan add 'smooth' and 'smearing_sigma' 2009-08-28
		ofs << 0 << " smooth" << std::endl;
		ofs << 0 << " smearing_sigma" << std::endl;

		ofs << Exx_Abfs::Jle::tolerence << " tolerence" << std::endl;

		ofs << Exx_Abfs::Jle::Lmax << " lmax" << std::endl;

		ofs << GlobalC::kv.nkstot << " nks" << std::endl;
		ofs	<< index_lcaos[TA].count_size * index_lcaos[TB].count_size << " nbands" << std::endl;
		
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
		assert( ecut_numberA == ecut_numberB );
		ofs	<< ecut_numberA << " ne" << std::endl;
		
		ofs << "<WEIGHT_OF_KPOINTS>" << std::endl;
		for( int ik=0; ik!=GlobalC::kv.nkstot; ++ik )		
		{
			ofs << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z;
			ofs << " " << GlobalC::kv.wk[ik] * 0.5 << std::endl;
		}
		ofs << "</WEIGHT_OF_KPOINTS>" << std::endl;

		ofs << std::endl;
	};
	
	auto print_Q = [&]( std::ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < jY | Psi >
		---------------------*/	

		ofs<< "<OVERLAP_Q>" << std::endl;
		
		for( int ik=0; ik!=GlobalC::kv.nkstot; ++ik )
			for( size_t LA=0; LA!=range_lcaos[TA].size(); ++LA )
				for( size_t NA=0; NA!=range_lcaos[TA][LA].N; ++NA )
					for( size_t MA=0; MA!=range_lcaos[TA][LA].M; ++MA )
						for( size_t LB=0; LB!=range_lcaos[TB].size(); ++LB )
							for( size_t NB=0; NB!=range_lcaos[TB][LB].N; ++NB )	
								for( size_t MB=0; MB!=range_lcaos[TB][LB].M; ++MB )
								{
									const std::vector<ModuleBase::matrix> & matrix_Q = matrixes_Q.at(TA).at(IA).at(TB).at(IB);
									const size_t index_lcao 
										= Exx_Abfs::Abfs_Index::get_index_index( 
											index_lcaos,TA,LA,NA,MA, 
											index_lcaos,TB,LB,NB,MB );
											
									auto f = [&]( size_t Tj, size_t Ij )
									{
										for( size_t Lj=0; Lj!=range_jles[Tj].size(); ++Lj )
											for( size_t Mj=0; Mj!=range_jles[Tj][Lj].M; ++Mj )				// attention! M first N second
												for( size_t Nj=0; Nj!=range_jles[Tj][Lj].N; ++Nj )
													ofs << matrix_Q[Ij]( 
															index_lcao, 
															index_jles[Tj][Lj][Nj][Mj] )
														 * scale << "\t" 
														<< 0.0 * scale << std::endl;
									};

									if( TA==TB && IA==IB )
									{
										f(TA,0);
									}
									else
									{
										f(TA,0);
										f(TB,1);
									}
								}

		
		ofs<< "</OVERLAP_Q>" << std::endl << std::endl;
	};


	auto print_S = [&]( std::ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < jY | jY >
		---------------------*/		
		auto q1q2 = [&]( size_t T1, size_t I1, size_t T2, size_t I2, size_t L1, size_t M1, size_t L2, size_t M2)
		{
			const ModuleBase::matrix & matrix_S = matrixes_S.at(T1).at(I1).at(T2).at(I2);
			for( size_t N1=0; N1!=range_jles[T1][L1].N; ++N1)
				for( size_t N2=0; N2!=range_jles[T2][L2].N; ++N2)
					ofs<< matrix_S( index_jles[T1][L1][N1][M1], index_jles[T2][L2][N2][M2] ) * scale << "\t" << 0.0 * scale << std::endl;
		};
		
		auto atom2 = [&]( size_t T1, size_t I1, size_t T2, size_t I2, size_t L1, size_t M1)
		{
			for( size_t L2=0; L2!=range_jles[T2].size(); ++L2 )	
				for( size_t M2=0; M2!=range_jles[T2][L2].M; ++M2 )
					q1q2(T1,I1,T2,I2,L1,M1,L2,M2);
		};
		
		auto atom1 = [&]( size_t T1, size_t I1)
		{
			for( size_t L1=0; L1!=range_jles[T1].size(); ++L1 )
				for( size_t M1=0; M1!=range_jles[T1][L1].M; ++M1 )
				{
					if( TA==TB && IA==IB )
					{
						atom2(T1,I1,TA,IA,L1,M1);
					}
					else
					{
						atom2(T1,I1,TA,IA,L1,M1);
						atom2(T1,I1,TB,IB,L1,M1);
					}
				}
		};	
		
		ofs<< "<OVERLAP_Sq>" <<std::endl;
		
		for( int ik=0; ik!=GlobalC::kv.nkstot; ++ik )
		{
			if( TA==TB && IA==IB )
			{
				atom1(TA,IA);
			}
			else
			{
				atom1(TA,IA);
				atom1(TB,IB);
			}
		}	
		
		ofs<< "</OVERLAP_Sq>" << std::endl << std::endl;
	};


	auto print_V = [&]( std::ofstream &ofs, const size_t TA, const size_t IA, const size_t TB, const size_t IB, const double scale=1 )
	{
		/*---------------------
		  < Psi | Psi >
		---------------------*/
		ofs << "<OVERLAP_V>" << std::endl;
		
		const ModuleBase::matrix & matrix_V = matrixes_V.at(TA).at(IA).at(TB).at(IB);
		
		for( int ik=0; ik!=GlobalC::kv.nkstot; ++ik )
	//		for ib = 0 to GlobalV::NBANDS
			for( size_t LA=0; LA!=range_lcaos[TA].size(); ++LA )
				for( size_t NA=0; NA!=range_lcaos[TA][LA].N; ++NA )
					for( size_t MA=0; MA!=range_lcaos[TA][LA].M; ++MA )
						for( size_t LB=0; LB!=range_lcaos[TB].size(); ++LB )
							for( size_t NB=0; NB!=range_lcaos[TB][LB].N; ++NB )
								for( size_t MB=0; MB!=range_lcaos[TB][LB].M; ++MB )
									ofs << matrix_V( index_lcaos[TA][LA][NA][MA], index_lcaos[TB][LB][NB][MB] ) * scale << std::endl;
		
		ofs << "</OVERLAP_V>" << std::endl << std::endl;
	};

	
	auto cal_R2 = []( const size_t TA, const size_t IA, const size_t TB, const size_t IB ) ->double	
	{
		Numerical_Orbital::set_position( GlobalC::ucell.atoms[TA].tau[IA], GlobalC::ucell.atoms[TB].tau[IB] );
		const double R = Numerical_Orbital::get_distance()*GlobalC::ucell.lat0;
		return R*R;
	};
	
	for( size_t TA=0; TA!=GlobalC::ucell.ntype; ++TA )
	{
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
		{
			for( size_t TB=TA; TB!=GlobalC::ucell.ntype; ++TB )
			{
				for( size_t IB=((TB==TA)?IA:0); IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					std::ofstream ofs(( file_name_prefix+"matrix_"+ModuleBase::GlobalFunc::TO_STRING(TA)+"_"+ModuleBase::GlobalFunc::TO_STRING(IA)+"_"+ModuleBase::GlobalFunc::TO_STRING(TB)+"_"+ModuleBase::GlobalFunc::TO_STRING(IB) ).c_str());
					print_header( ofs, TA, IA, TB, IB );
//					const double scale = 1.0 / max( matrixes_V.at(TA).at(IA).at(TB).at(IB) );
					const double scale = 1;		// Peize Lin test
//					const double scale = cal_R2(TA,IA,TB,IB);
					print_Q( ofs, TA, IA, TB, IB, sqrt(scale) );
					print_S( ofs, TA, IA, TB, IB );
					print_V( ofs, TA, IA, TB, IB, scale );
				}
			}
		}
	}
}
