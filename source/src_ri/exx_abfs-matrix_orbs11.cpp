#include "exx_abfs-matrix_orbs11.h"

#include <algorithm>
#include <set>
#include "../src_pw/global.h"
#include "../module_orbital/ORB_read.h"
#include "../module_base/ylm.h"
#include "../module_base/global_function.h"

#include<sys/time.h>					// Peize Lin test
#include "../src_external/src_test/test_function.h"			// Peize Lin test 2016-04-05
#include "../src_external/src_test/src_ri/exx_lcao-test.h"
#include "../src_lcao/global_fp.h"

void Exx_Abfs::Matrix_Orbs11::init(
	const int mode, 
	const double kmesh_times, 
	const double rmesh_times)
{
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init");
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
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init::MOT.allocate\t"<<time_during(t_start)<<std::endl;
	int Lmax_used, Lmax;
//gettimeofday( &t_start, NULL);
	MOT.init_Table_Spherical_Bessel (2, mode, Lmax_used, Lmax, Exx_Abfs::Lmax,GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
//	MOT.init_OV_Tpair();							// for MOT.OV_L2plus1
//	MOT.Destroy_Table_Spherical_Bessel (Lmax_used);				// why?
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init::MOT.init_Table_Spherical_Bessel\t"<<time_during(t_start)<<std::endl;

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	//liaochen add 2010/4/29
//gettimeofday( &t_start, NULL);
	ModuleBase::Ylm::set_coefficients ();
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init::Ylm\t"<<time_during(t_start)<<std::endl;
	
	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
//gettimeofday( &t_start, NULL);
	MGT.init_Gaunt_CH( Lmax );
	MGT.init_Gaunt( Lmax );
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init::MGT\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}

void Exx_Abfs::Matrix_Orbs11::init_radial(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A, 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B)
{ 
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init_radial");
	for( size_t TA = 0; TA!=orb_A.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.size(); ++TB )
			for( int LA=0; LA!=orb_A[TA].size(); ++LA )
				for( size_t NA=0; NA!=orb_A[TA][LA].size(); ++NA )
					for( int LB=0; LB!=orb_B[TB].size(); ++LB )
						for( size_t NB=0; NB!=orb_B[TB][LB].size(); ++NB )
							center2_orb11_s[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb11(
									orb_A[TA][LA][NA], 
									orb_B[TB][LB][NB],
									MOT, MGT)));
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init_radial\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}

void Exx_Abfs::Matrix_Orbs11::init_radial(
	const LCAO_Orbitals &orb_A, 
	const LCAO_Orbitals &orb_B)
{ 
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init_radial");
	for( size_t TA = 0; TA!=orb_A.get_ntype(); ++TA )
		for( size_t TB=0; TB!=orb_B.get_ntype(); ++TB )
			for( int LA=0; LA<=orb_A.Phi[TA].getLmax(); ++LA )
				for( size_t NA=0; NA!=orb_A.Phi[TA].getNchi(LA); ++NA )
					for( int LB=0; LB<=orb_B.Phi[TB].getLmax(); ++LB )
						for( size_t NB=0; NB!=orb_B.Phi[TB].getNchi(LB); ++NB )
							center2_orb11_s[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb11(
									orb_A.Phi[TA].PhiLN(LA,NA),								
									orb_B.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init_radial\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}

void Exx_Abfs::Matrix_Orbs11::init_radial_table()
{
//std::ofstream ofs(exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
//timeval t_start;
//gettimeofday( &t_start, NULL);
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init_radial_table");
	for( auto &coA : center2_orb11_s )
		for( auto &coB : coA.second )
			for( auto &coC : coB.second )
				for( auto &coD : coC.second )
					for( auto &coE : coD.second )
						for( auto &coF : coE.second )
							coF.second.init_radial_table();
//ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init_radial_table\t"<<time_during(t_start)<<std::endl;
//ofs.close();
}

void Exx_Abfs::Matrix_Orbs11::init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs )
{	
std::ofstream ofs(GlobalC::exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init_radial_table_Rs");
	for( const auto &RsA : Rs )
		for( const auto &RsB : RsA.second )
		{
			if( auto* const center2_orb11_sAB = static_cast<std::map<int,std::map<size_t,std::map<int,std::map<size_t,Center2_Orb::Orb11>>>>*const>(
						ModuleBase::GlobalFunc::MAP_EXIST(center2_orb11_s, RsA.first, RsB.first)) )
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
				for( auto &coC : *center2_orb11_sAB )
					for( auto &coD : coC.second )
						for( auto &coE : coD.second )
							for( auto &coF : coE.second )
								coF.second.init_radial_table(radials);
ofs<<time_during(t_small)<<std::endl;
			}
		}
ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::init_radial_table_Rs\t"<<time_during(t_start)<<std::endl;
ofs.close();
}

/*
void Exx_Abfs::Matrix_Orbs11::init_radial_table()
{
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","init_radial_table");
	for( auto &co1 : center2_orb11_s )
	{
		const size_t TA = co1.first;
		for( auto &co2 : co1.second )
		{
			const size_t TB = co2.first;
			for( auto &co3 : co2.second )
			{
				const size_t LA = co3.first;
				for( auto &co4 : co3.second )
				{
					const size_t NA = co4.first;
					for( auto &co5 : co4.second )
					{
						const size_t LB = co5.first;
						for( auto &co6 : co5.second )
						{
							const size_t NB = co6.first;
//if(TA==TB&&LA==LB&&NA==NB)	{ exx_cout_flag=true; std::cout<<TA<<"\t"<<LA<<"\t"<<NA<<std::endl; }
							co6.second.init_radial_table();
//exx_cout_flag=false;
						}
					}
				}
			}
		}
	}
}
*/
/*
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> Matrix_Orbs11::cal_overlap_matrix( 
	const Exx_Abfs::Abfs_Index::Index &index_r, 
	const Exx_Abfs::Abfs_Index::Index &index_c )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> matrix_V;
	
	for( auto &co1 : center2_orb11_s )
	{
		const size_t TA = co1.first;
		for (size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );
			GlobalC::GridD.Find_atom(tauA);

			for( auto &co2 : co1.second )
			{
				const int LA = co2.first;
				for( size_t MA=0; MA!=2*LA+1; ++MA)
				{
					for( auto &co3 : co2.second )
					{
						const size_t NA = co3.first;
						for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
						{
							for( auto &co4 : co3.second )
							{
								const size_t TB = co4.first;
								if( TB != GlobalC::GridD.getType(ad) )
									continue;
								const ModuleBase::Vector3<double> &tauB( GlobalC::GridD.getAdjacentTau(ad) );
								const size_t IB = GlobalC::GridD.getNatom(ad);
																					
								for( auto &co5 : co4.second )
								{
									const int LB = co5.first;
									for( size_t MB=0; MB!=2*LB+1; ++MB)
									{
										for( auto &co6 : co5.second )
										{
											const size_t NB = co6.firs;
											co6.second.cal_ST_Phi12_R(tauA, tauB, MA, MB);
											
											matrix_V[TA][IA][TB][IB]( index_r[TA][LA][MA][NA], index_c[TB][LB][MB][NB] ) = co6.second.olm[0];
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
	return matrix_V;
}
*/

ModuleBase::matrix Exx_Abfs::Matrix_Orbs11::cal_overlap_matrix( 
	const size_t TA, 
	const size_t TB, 
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_r, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_c ) const
{
//	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","cal_overlap_matrix");
	
	ModuleBase::matrix m( index_r[TA].count_size, index_c[TB].count_size );
	
	for( const auto &co3 : center2_orb11_s.at(TA).at(TB) )
	{
		const int LA = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA = co4.first;
			for( size_t MA=0; MA!=2*LA+1; ++MA )
			{
				for( const auto &co5 : co4.second )
				{
					const int LB = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NB = co6.first;	
						for( size_t MB=0; MB!=2*LB+1; ++MB )
						{
//if(TA==TB&&LA==LB&&NA==NB&&MA==MB)	{exx_cout_flag=true; std::cout<<TA<<" "<<LA<<" "<<NA<<" "<<MA<<" "<<tauA<<" "<<tauB<<std::endl; }
							m( index_r[TA][LA][NA][MA], index_c[TB][LB][NB][MB] ) 
							= co6.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA, MB );
//exx_cout_flag=false;
						}
					}
				}
			}
		}
	}
	return m;
}

std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> Exx_Abfs::Matrix_Orbs11::cal_overlap_matrix( 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_r, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_c ) const
{
std::ofstream ofs(GlobalC::exx_lcao.test_dir.process+"time_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK),std::ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	ModuleBase::TITLE("Exx_Abfs::Matrix_Orbs11","cal_overlap_matrix");
	
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> matrixes;
	
	for( const auto &co1 : center2_orb11_s )
	{
		const size_t TA = co1.first;
		for (size_t IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );

			for( const auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for (size_t IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB)
				{
					const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[TB].tau[IB] );

					matrixes[TA][IA][TB][IB] = cal_overlap_matrix( TA, TB, tauA, tauB, index_r, index_c );
				}
			}
		}
	}
ofs<<"TIME@Exx_Abfs::Matrix_Orbs11::cal_overlap_matrix\t"<<time_during(t_start)<<std::endl;
ofs.close();
	return matrixes;
}
