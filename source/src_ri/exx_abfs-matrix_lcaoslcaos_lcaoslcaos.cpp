#include "exx_abfs-matrix_lcaoslcaos_lcaoslcaos.h"
#include "../src_pw/global.h"
#include "../module_orbital/ORB_read.h"
#include "../module_base/ylm.h"
#include "../src_external/src_test/test_function.h"			// Peize Lin test 2016-04-05

void Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::init(
	const int mode, 
	const double kmesh_times, 
	const double rmesh_times)
{
	//=========================================
	// (1) MOT: make overlap table.
	//=========================================
	MOT.allocate(
		GlobalC::ORB.get_ntype(),							// number of atom types
		GlobalC::ORB.get_lmax(),								// max L used to calculate overlap
		static_cast<int>(GlobalC::ORB.get_kmesh() * kmesh_times) | 1,			// kpoints, for integration in k space
		GlobalC::ORB.get_Rmax() * rmesh_times,				// max value of radial table
		GlobalC::ORB.get_dR(),								// delta R, for making radial table
//		GlobalC::ORB.get_dk() / kmesh_times);				// delta k, for integration in k space
		GlobalC::ORB.get_dk());											// Peize Lin change 2017-04-16
	int Lmax_used, Lmax;
	MOT.init_Table_Spherical_Bessel (4,mode, Lmax_used, Lmax, Exx_Abfs::Lmax,GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
//	MOT.init_OV_Tpair();							// for MOT.OV_L2plus1
//	MOT.Destroy_Table_Spherical_Bessel (Lmax_used);				// why?

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	ModuleBase::Ylm::set_coefficients ();
	
	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
	MGT.init_Gaunt_CH( 2*Lmax+1 );			// why +1
	MGT.init_Gaunt( 2*Lmax+1 );
}

void Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::init_radial( 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A, 
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B )
{
	for( size_t TA = 0; TA!=orb_A.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.size(); ++TB )
			for( size_t LA=0; LA!=orb_A[TA].size(); ++LA )
				for( size_t NA=0; NA!=orb_A[TA][LA].size(); ++NA )
					for( size_t LB=0; LB!=orb_B[TB].size(); ++LB )
						for( size_t NB=0; NB!=orb_B[TB][LB].size(); ++NB )
							center2_orb22_s[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb22(
									orb_A[TA][LA][NA],
									orb_A[TA][LA][NA],									
									orb_B[TB][LB][NB],
									orb_B[TB][LB][NB],
									MOT, MGT)));
}

void Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::init_radial( 
	const LCAO_Orbitals &orb_A, 
	const LCAO_Orbitals &orb_B )
{
	for( size_t TA = 0; TA!=orb_A.get_ntype(); ++TA )
		for( size_t TB=0; TB!=orb_B.get_ntype(); ++TB )
			for( size_t LA=0; LA<=orb_A.Phi[TA].getLmax(); ++LA )
				for( size_t NA=0; NA!=orb_A.Phi[TA].getNchi(LA); ++NA )
					for( size_t LB=0; LB<=orb_B.Phi[TB].getLmax(); ++LB )
						for( size_t NB=0; NB!=orb_B.Phi[TB].getNchi(LB); ++NB )
							center2_orb22_s[TA][TB][LA][NA][LB].insert( 
								make_pair(NB, Center2_Orb::Orb22(
									orb_A.Phi[TA].PhiLN(LA,NA),
									orb_A.Phi[TA].PhiLN(LA,NA),									
									orb_B.Phi[TB].PhiLN(LB,NB),
									orb_B.Phi[TB].PhiLN(LB,NB),
									MOT, MGT)));
}

void Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::init_radial_table()
{
	for( auto &co1 : center2_orb22_s )
		for( auto &co2 : co1.second )
			for( auto &co3 : co2.second )
				for( auto &co4 : co3.second )
					for( auto &co5 : co4.second )
						for( auto &co6 : co5.second )
							co6.second.init_radial_table();
}

void Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::init_radial_table( std::map<size_t,std::map<size_t,set<double>>> &Rs )
{
	for( auto &co1 : center2_orb22_s )
		for( auto &co2 : co1.second )
		{
			set<size_t> radials;
			for( const double &R : Rs[co1.first][co2.first] )
			{
				const double position = R * GlobalC::ucell.lat0 / MOT.dr;
				const size_t iq = static_cast<size_t>(position);
				for( size_t i=0; i!=4; ++i )
					radials.insert(iq+i);
			}
			
			for( auto &co3 : co2.second )
				for( auto &co4 : co3.second )
					for( auto &co5 : co4.second )
						for( auto &co6 : co5.second )
							co6.second.init_radial_table(radials);
		}
}

/*void Matrix_Phiphi_Phiphi::cal_overlap_matrix()
{
	for( auto &co1 : center2_orb22_s )
	{
		const size_t T = co1.first;
		for (size_t IA=0; IA!=GlobalC::ucell.atoms[T].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[T].tau[IA] );
			GlobalC::GridD.Find_atom(tauA);

			for( auto &co2 : co1.second )
			{
				const size_t LA = co2.first;
				for( size_t MA=0; MA!=2*LA+1; ++MA)
				{
					for( auto &co3 : co2.second )
					{

						for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
						{
							if( T != GlobalC::GridD.getType(ad) )
								continue;
							const ModuleBase::Vector3<double> &tauB( GlobalC::GridD.getAdjacentTau(ad) );

							for( auto &co4 : co3.second )
							{
								const size_t LB = co4.first;
								for( size_t MB=0; MB!=2*LB+1; ++MB)
								{
									for( auto &co5 : co4.second )
									{
										co5.second.cal_ST_Phi12_R(tauA, tauB, MA, MA, MB, MB);

									}
								}

							}
						}
					}
				}
			}
		}
	}
}*/

std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> Exx_Abfs::Matrix_Lcaoslcaos_Lcaoslcaos::cal_overlap_matrix(
	const ModuleBase::Element_Basis_Index::IndexLNM &index_r, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_c )
{
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,ModuleBase::matrix>>>> matrixes;

	for( auto &co1 : center2_orb22_s )
	{
		const int TA = co1.first;
		for (int IA=0; IA!=GlobalC::ucell.atoms[TA].na; ++IA)
		{
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[TA].tau[IA] );

			for( auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for ( int IB=0; IB!=GlobalC::ucell.atoms[TB].na; ++IB )
				{
					const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[TB].tau[IB] );	
					
					matrixes[TA][IA][TB][IB].create( index_r[TA].count_size, index_c[TB].count_size );			
														
					for( auto &co3 : co2.second )
					{
						const int LA = co3.first;
						for( auto &co4 : co3.second )
						{
							const size_t NA = co4.first;
							for( int MA=0; MA!=2*LA+1; ++MA)
							{
								for( auto &co5 : co4.second )
								{
									const int LB = co5.first;
									for( auto &co6 : co5.second )
									{
										const size_t NB = co6.first;
										for( int MB=0; MB!=2*LB+1; ++MB)
										{
											matrixes[TA][IA][TB][IB]( index_r[TA][LA][NA][MA], index_c[TB][LB][NB][MB] ) 
											= co6.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA, MA, MB, MB );
											
											// Peize Lin test
//											std::cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<"\t"<<LA<<"\t"<<MA<<"\t"<<NA<<"\t"<<LB<<"\t"<<MB<<"\t"<<NB<<"\t"<<std::endl; //co6.second.olm[0]<<std::endl;
										}
									}
								}
							}
						}
					}
					// Peize Lin test
//					std::cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<std::endl;
//					std::cout<<matrixes[TA][IA][TB][IB]<<std::endl;
				}				
			}
		}
	}
	return matrixes;
}

/*void Matrix_Phiphi_Phiphi::cal_overlap_matrix()
{
	for( auto &co1 : center2_orb22_s )
	{
		const int T = co1.first;

// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<T<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}		
		for (int IA=0; IA!=GlobalC::ucell.atoms[T].na; ++IA)
		{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<" "<<IA<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}			
			const ModuleBase::Vector3<double> &tauA( GlobalC::ucell.atoms[T].tau[IA] );

			for( auto &co2 : co1.second )
			{
				const int LA = co2.first;
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"  "<<LA<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}				
				for( int MA=0; MA!=2*LA+1; ++MA)
				{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"   "<<MA<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}					
					for( auto &co3 : co2.second )
					{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"    "<<co3.first<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}						
						for (int IB=0; IB!=GlobalC::ucell.atoms[T].na; ++IB)
						{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"     "<<IB<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}							
							const ModuleBase::Vector3<double> &tauB( GlobalC::ucell.atoms[T].tau[IB] );

							for( auto &co4 : co3.second )
							{
								const int LB = co4.first;
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"      "<<LB<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}								
								for( int MB=0; MB!=2*LB+1; ++MB)
								{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"       "<<MB<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}								
									for( auto &co5 : co4.second )
									{
// Peize Lin test 2015-04-05
{
	std::stringstream ss;
	ss<<"        "<<co5.first<<std::endl;
	MPI_RANK_OFSTREAM( "Matrix_Phiphi_Phiphi::cal_overlap_matrix", ss);
}								
										co5.second.cal_ST_Phi12_R(tauA, tauB, MA, MA, MB, MB);
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
*/
