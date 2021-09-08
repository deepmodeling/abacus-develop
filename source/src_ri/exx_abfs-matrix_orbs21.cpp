#include "exx_abfs-matrix_orbs21.h"

#include "exx_abfs.h"
#include "exx_abfs-abfs_index.h"

#include <algorithm>
#include <set>
#include "../src_pw/global.h"
#include "../module_orbital/ORB_read.h"
#include "../module_base/ylm.h"

#include<sys/time.h>					// Peize Lin test
#include "../src_external/src_test/test_function.h"			// Peize Lin test 2016-04-05
#include "../src_external/src_test/src_ri/exx_lcao-test.h"
#include "../src_lcao/global_fp.h"

void Exx_Abfs::Matrix_Orbs21::init(
	const int mode,
	const double kmesh_times,
	const double rmesh_times)
{
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","init");
	//=========================================
	// (1) MOT: make overlap table.
	//=========================================
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);
	MOT.allocate(
		GlobalC::ORB.get_ntype(),							// number of atom types
		std::max( GlobalC::ORB.get_lmax(), Exx_Abfs::Lmax ),	// max L used to calculate overlap
		static_cast<int>(GlobalC::ORB.get_kmesh() * kmesh_times) | 1,				// kpoints, for integration in k space
		GlobalC::ORB.get_Rmax() * rmesh_times,				// max value of radial table
		GlobalC::ORB.get_dR(),								// delta R, for making radial table
//		GlobalC::ORB.get_dk() / kmesh_times);				// delta k, for integration in k space
		GlobalC::ORB.get_dk());											// Peize Lin change 2017-04-16
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init::MOT.allocate\t"<<time_during(t_start)<<std::endl;
	int Lmax_used, Lmax;
//gettimeofday( &t_start, NULL);
	MOT.init_Table_Spherical_Bessel (3,mode, Lmax_used, Lmax, Exx_Abfs::Lmax,GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
//	MOT.init_OV_Tpair();							// for MOT.OV_L2plus1
//	MOT.Destroy_Table_Spherical_Bessel (Lmax_used);				// why?
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init::MOT.init_Table_Spherical_Bessel\t"<<time_during(t_start)<<std::endl;

	//=========================================
	// (2) init Ylm Coef
	//=========================================
//gettimeofday( &t_start, NULL);
	ModuleBase::Ylm::set_coefficients ();
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init::Ylm\t"<<time_during(t_start)<<std::endl;

	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
//gettimeofday( &t_start, NULL);
	MGT.init_Gaunt_CH( 2*Lmax+1 );			// why +1
	MGT.init_Gaunt( 2*Lmax+1 );
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init::MGT\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}



void Exx_Abfs::Matrix_Orbs21::init_radial(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A2,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B )
{
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);

	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","init_radial");

	assert(orb_A1.size()==orb_A2.size());
	for( size_t TA=0;  TA!=orb_A1.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.size(); ++TB )
			for( int LA1=0; LA1!=orb_A1[TA].size(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1[TA][LA1].size(); ++NA1 )
					for( int LA2=0; LA2!=orb_A2[TA].size(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2[TA][LA2].size(); ++NA2 )
							for( int LB=0; LB!=orb_B[TB].size(); ++LB )
								for( size_t NB=0; NB!=orb_B[TB][LB].size(); ++NB )
									center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
										make_pair(NB, Center2_Orb::Orb21(
											orb_A1[TA][LA1][NA1],
											orb_A2[TA][LA2][NA2],
											orb_B[TB][LB][NB],
											MOT, MGT)));
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init_radial\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}


void Exx_Abfs::Matrix_Orbs21::init_radial(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1,
	const LCAO_Orbitals &orb_A2,
	const LCAO_Orbitals &orb_B )
{
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);

	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","init_radial");

	assert( orb_A1.size() == orb_A2.get_ntype() );
	for( size_t TA=0;  TA!=orb_A1.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.get_ntype(); ++TB)
			for( int LA1=0; LA1!=orb_A1[TA].size(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1[TA][LA1].size(); ++NA1 )
					for( int LA2=0; LA2<=orb_A2.Phi[TA].getLmax(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2.Phi[TA].getNchi(LA2); ++NA2 )
							for( int LB=0; LB<=orb_B.Phi[TB].getLmax(); ++LB )
								for( size_t NB=0; NB!=orb_B.Phi[TB].getNchi(LB); ++NB )
									center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
										make_pair(NB, Center2_Orb::Orb21(
											orb_A1[TA][LA1][NA1],
											orb_A2.Phi[TA].PhiLN(LA2,NA2),
											orb_B.Phi[TB].PhiLN(LB,NB),
											MOT, MGT)));
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init_radial\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}


void Exx_Abfs::Matrix_Orbs21::init_radial_table()
{
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);

	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","init_radial_table");

	for( auto &coA : center2_orb21_s )
		for( auto &coB : coA.second )
			for( auto &coC : coB.second )
				for( auto &coD : coC.second )
					for( auto &coE : coD.second )
						for( auto &coF : coE.second )
							for( auto &coG : coF.second )
								for( auto &coH : coG.second )
									coH.second.init_radial_table();
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init_radial_table\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}

void Exx_Abfs::Matrix_Orbs21::init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs )
{
std::ofstream ofs(GlobalC::exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;
gettimeofday(&t_start, NULL);

	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","init_radial_table_Rs");

	for( const auto &RsA : Rs )
		for( const auto &RsB : RsA.second )
		{
			if( auto* const center2_orb21_sAB = static_cast<std::map<int,std::map<size_t,std::map<int,std::map<size_t,std::map<int,std::map<size_t,Center2_Orb::Orb21>>>>>>*const>(
						ModuleBase::GlobalFunc::MAP_EXIST(center2_orb21_s, RsA.first, RsB.first)) )
			{
timeval t_small;
gettimeofday(&t_small, NULL);
				std::set<size_t> radials;
				for( const double &R : RsB.second )
				{
					const double position = R * GlobalC::ucell.lat0 / MOT.dr;
					const size_t iq = static_cast<size_t>(position);
					for( size_t i=0; i!=4; ++i )
						radials.insert(iq+i);
				}
ofs<<"\t"<<RsA.first<<"\t"<<RsB.first<<"\t"<<time_during(t_small)<<"\t"<<flush;
gettimeofday(&t_small, NULL);
				for( auto &coC : *center2_orb21_sAB )
					for( auto &coD : coC.second )
						for( auto &coE : coD.second )
							for( auto &coF : coE.second )
								for( auto &coG : coF.second )
									for( auto &coH : coG.second )
										coH.second.init_radial_table(radials);
ofs<<time_during(t_small)<<std::endl;
			}
		}
ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::init_radial_table_Rs\t"<<time_during(t_start)<<std::endl;
ofs.close();
}

/*
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Matrix_Abfsphi_Phi::cal_overlap_matrix(
	const Exx_Abfs::Abfs_Index::Index &index_abfs,
	const Exx_Abfs::Abfs_Index::Index &index_orb )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> matrix_A;

	for( auto &co1 : center2_orb21_s )
	{
		const size_t TA = co1.first;
		for (size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );
			GlobalC::GridD.Find_atom(tauA);

			for( auto &co2 : co1.second )
			{
				const int LA1 = co2.first;
				for( size_t MA1=0; MA1!=2*LA1+1; ++MA1)
				{
					for( auto &co3 : co2.second )
					{
						const size_t NA1 = co3.first;
						for( auto &co4 : co3.second )
						{
							const int LA2 = co4.first;
							for( size_t MA2=0; MA2!=2*LA2+1; ++MA2)
							{
								for( auto &co5 : co4.second )
								{
									const size_t NA2 = co4.first;
									for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
									{
										for( auto &co6 : co5.second )
										{
											const size_t TB = co6.first;
											if( TB != GlobalC::GridD.getType(ad) )
												continue;
											const ModuleBase::Vector3<double> &tauB( GlobalC::GridD.getAdjacentTau(ad) );
											const size_t IB = GlobalC::GridD.getNatom(ad);

											for( auto &co7 : co6.second )
											{
												const int LB = co7.first;
												for( size_t MB=0; MB!=2*LB+1; ++MB)
												{
													for( auto &co8 : co7.second )
													{
														const size_t NB = co8.first;
														co8.second.cal_ST_Phi12_R(tauA, tauB, MA1, MA2, MB);

														matrix_A[TA][IA][TB][IB][0]( Exx_Abfs::Abfs_Index::get_index_index( index_orb,TA,LA2,MA2,NA2, index_orb,TB,LB,MB,NB ), index_abfs[TA][LA1][MA1][NA1] )
														= matrix_A[TB][IB][TA][IA][1]( Exx_Abfs::Abfs_Index::get_index_index( index_orb,TB,LB,MB,NB, index_orb,TA,LA2,MA2,NA2 ), index_abfs[TA][LA1][MA1][NA1] )
														= co8.second.olm[0];
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return matrix_A;
}*/

ModuleBase::matrix Exx_Abfs::Matrix_Orbs21::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order) const
{
	ModuleBase::matrix m;
	switch(matrix_order)
	{
		case Matrix_Order::A2B_A1:	m.create( index_A2[TA].count_size*index_B [TB].count_size, index_A1[TA].count_size );	break;
		case Matrix_Order::BA2_A1:	m.create( index_B [TB].count_size*index_A2[TA].count_size, index_A1[TA].count_size );	break;
		default:	throw std::invalid_argument( "Matrix_Order wrong in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__) );
	}

	for( const auto &co3 : center2_orb21_s.at(TA).at(TB) )
	{
		const int LA1 = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA1 = co4.first;
			for( size_t MA1=0; MA1!=2*LA1+1; ++MA1 )
			{
				for( const auto &co5 : co4.second )
				{
					const int LA2 = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NA2 = co6.first;
						for( size_t MA2=0; MA2!=2*LA2+1; ++MA2 )
						{
							for( const auto &co7 : co6.second )
							{
								const int LB = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB = co8.first;
									for( size_t MB=0; MB!=2*LB+1; ++MB )
									{
										const double overlap = co8.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB );

										switch(matrix_order)
										{
											case Matrix_Order::A2B_A1:
												m(
													Exx_Abfs::Abfs_Index::get_index_index(
														index_A2,TA,LA2,NA2,MA2,
														index_B,TB,LB,NB,MB ),
													index_A1[TA][LA1][NA1][MA1] )
												= overlap;
												break;
											case Matrix_Order::BA2_A1:
												m(
													Exx_Abfs::Abfs_Index::get_index_index(
														index_B,TB,LB,NB,MB,
														index_A2,TA,LA2,NA2,MA2),
													index_A1[TA][LA1][NA1][MA1] )
												= overlap;
												break;
											default:	throw std::invalid_argument( "Matrix_Order wrong in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__) );
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return m;
}

std::vector<ModuleBase::matrix> Exx_Abfs::Matrix_Orbs21::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B) const
{
	ModuleBase::matrix m_A2B_A1( index_A2[TA].count_size*index_B[TB].count_size, index_A1[TA].count_size );
	ModuleBase::matrix m_BA2_A1( index_B[TB].count_size*index_A2[TA].count_size, index_A1[TA].count_size );

	for( const auto &co3 : center2_orb21_s.at(TA).at(TB) )
	{
		const int LA1 = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA1 = co4.first;
			for( size_t MA1=0; MA1!=2*LA1+1; ++MA1 )
			{
				for( const auto &co5 : co4.second )
				{
					const int LA2 = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NA2 = co6.first;
						for( size_t MA2=0; MA2!=2*LA2+1; ++MA2 )
						{
							for( const auto &co7 : co6.second )
							{
								const int LB = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB = co8.first;
									for( size_t MB=0; MB!=2*LB+1; ++MB )
									{
										const double overlap = co8.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA1, MA2, MB );

										m_A2B_A1(
											Exx_Abfs::Abfs_Index::get_index_index(
												index_A2,TA,LA2,NA2,MA2,
												index_B,TB,LB,NB,MB ),
											index_A1[TA][LA1][NA1][MA1] )
										= m_BA2_A1(
											Exx_Abfs::Abfs_Index::get_index_index(
												index_B,TB,LB,NB,MB,
												index_A2,TA,LA2,NA2,MA2),
											index_A1[TA][LA1][NA1][MA1] )
										= overlap;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return std::vector<ModuleBase::matrix>{ std::move(m_A2B_A1), std::move(m_BA2_A1) };
}

std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> Exx_Abfs::Matrix_Orbs21::cal_overlap_matrix(
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A1,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A2,
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B) const
{
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);

	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs21","cal_overlap_matrix");

	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,std::vector<ModuleBase::matrix>>>>> matrixes;

	for( const auto &co1 : center2_orb21_s )
	{
		const size_t TA = co1.first;
		for( size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA )
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );

			for( const auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for( size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[TB].tau[IB] );

					const std::vector<ModuleBase::matrix> &&m = cal_overlap_matrix( TA, TB, tauA, tauB, index_A1, index_A2, index_B );
					matrixes[TA][IA][TB][IB].resize(2);
					matrixes[TA][IA][TB][IB][0] = std::move(m[0]);
					matrixes[TB][IB][TA][IA].resize(2);
					matrixes[TB][IB][TA][IA][1] = std::move(m[1]);
/*
					matrixes[TA][IA][TB][IB].resize(2);
					matrixes[TA][IA][TB][IB][0].create( index_A2[TA].count_size*index_B[TB].count_size, index_A1[TA].count_size );
					matrixes[TB][IB][TA][IA].resize(2);
					matrixes[TB][IB][TA][IA][1].create( index_B[TB].count_size*index_A2[TA].count_size, index_A1[TA].count_size );

					for( const auto &co3 : co2.second )
					{
						const size_t LA1 = co3.first;
						for( const auto &co4 : co3.second )
						{
							const size_t NA1 = co4.first;
							for( size_t MA1=0; MA1!=2*LA1+1; ++MA1 )
							{
								for( const auto &co5 : co4.second )
								{
									const size_t LA2 = co5.first;
									for( const auto &co6 : co5.second )
									{
										const size_t NA2 = co6.first;
										for( size_t MA2=0; MA2!=2*LA2+1; ++MA2 )
										{
											for( const auto &co7 : co6.second )
											{
												const size_t LB = co7.first;
												for( const auto &co8 : co7.second )
												{
													const size_t NB = co8.first;
													for( size_t MB=0; MB!=2*LB+1; ++MB )
													{

														// Peize Lin test											// Peize Lin test
//														{
//															std::ofstream ofs("orb11-i_exp",std::ofstream::app);
//															ofs<<IA<<"\t"<<IB<<"\t"<<LA1<<"\t"<<MA1<<"\t"<<NA1<<"\t"<<LA2<<"\t"<<MA2<<"\t"<<NA2<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<std::endl;
//															ofs.close();
//														}
//														{
//															std::ofstream ofs("orb11-Gaunt_solid_A_B_AB",std::ofstream::app);
//															ofs<<IA<<"\t"<<IB<<"\t"<<LA1<<"\t"<<MA1<<"\t"<<NA1<<"\t"<<LA2<<"\t"<<MA2<<"\t"<<NA2<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<std::endl;
//															ofs.close();
//														}
//														{
//															std::ofstream ofs("orb11-Gaunt_A_B_AB",std::ofstream::app);
//															ofs<<IA<<"\t"<<IB<<"\t"<<LA1<<"\t"<<MA1<<"\t"<<NA1<<"\t"<<LA2<<"\t"<<MA2<<"\t"<<NA2<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<std::endl;
//															ofs.close();
//														}
//														{
//															std::ofstream ofs("orb11-Interp_Tlm",std::ofstream::app);
//															ofs<<IA<<"\t"<<IB<<"\t"<<LA1<<"\t"<<MA1<<"\t"<<NA1<<"\t"<<LA2<<"\t"<<MA2<<"\t"<<NA2<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<std::endl;
//															ofs.close();
//														}
//														{
//															std::ofstream ofs("orb11-rly",std::ofstream::app);
//															ofs<<IA<<"\t"<<IB<<"\t"<<LA1<<"\t"<<MA1<<"\t"<<NA1<<"\t"<<LA2<<"\t"<<MA2<<"\t"<<NA2<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<std::endl;
//															ofs.close();
//														}

														matrixes[TA][IA][TB][IB][0](
															Exx_Abfs::Abfs_Index::get_index_index(
																index_A2,TA,LA2,NA2,MA2,
																index_B,TB,LB,NB,MB ),
															index_A1[TA][LA1][NA1][MA1] )
														= matrixes[TB][IB][TA][IA][1](
															Exx_Abfs::Abfs_Index::get_index_index(
																index_B,TB,LB,NB,MB,
																index_A2,TA,LA2,NA2,MA2 ),
															index_A1[TA][LA1][NA1][MA1] )
														= co8.second.cal_overlap( tauA, tauB, MA1, MA2, MB );
													}
												}
											}
										}
									}
								}
							}
						}
					}
*/
				}
			}
		}
	}

	// matrixes[T][I][T][I][0] = matrixes[T][I][T][I][1], so delete repeat
	for( auto m1 : matrixes )
	{
		const size_t T = m1.first;
		for( auto m2 : m1.second )
		{
			const size_t I = m2.first;
			matrixes[T][I][T][I].resize(1);
		}
	}
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs21::cal_overlap_matrix\t"<<time_during(t_start)<<std::endl;
//ofs.close();

	return matrixes;
}
