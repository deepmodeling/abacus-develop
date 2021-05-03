#include "exx_abfs-matrix_orbs22.h"

#include "exx_abfs.h"
#include "exx_abfs-abfs_index.h"

#include "src_pw/global.h"
#include "src_lcao/ORB_read.h"
#include "src_global/ylm.h"

#include "src_external/src_test/test_function.h"			// Peize Lin test 2016-04-05
#include "src_external/src_test/src_ri/exx_lcao-test.h"
#include "src_lcao/global_fp.h"

void Exx_Abfs::Matrix_Orbs22::init(
	const int mode,
	const double kmesh_times,
	const double rmesh_times)
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	
	//=========================================
	// (1) MOT: make overlap table.
	//=========================================
	MOT.allocate(
		ORB.get_ntype(),							// number of atom types
		ORB.get_lmax(),								// max L used to calculate overlap
		static_cast<int>(ORB.get_kmesh() * kmesh_times) | 1,			// kpoints, for integration in k space
		ORB.get_Rmax() * rmesh_times,				// max value of radial table
		ORB.get_dR(),								// delta R, for making radial table
//		ORB.get_dk() / kmesh_times);				// delta k, for integration in k space
		ORB.get_dk());											// Peize Lin change 2017-04-16
	int Lmax_used, Lmax;
	MOT.init_Table_Spherical_Bessel (4,mode, Lmax_used, Lmax, Exx_Abfs::Lmax);
//	MOT.init_OV_Tpair();							// for MOT.OV_L2plus1
//	MOT.Destroy_Table_Spherical_Bessel (Lmax_used);				// why?

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	Ylm::set_coefficients ();

	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
	MGT.init_Gaunt_CH( 2*Lmax+1 );			// why +1
	MGT.init_Gaunt( 2*Lmax+1 );
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::init\t"<<time_during(t_start)<<endl;
ofs.close();
}

void Exx_Abfs::Matrix_Orbs22::init_radial(
	const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A1,
	const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_A2,
	const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B1,
	const vector<vector<vector<Numerical_Orbital_Lm>>> &orb_B2 )
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	TITLE("Exx_Abfs::Matrix_Orbs22","init_radial");

	assert(orb_A1.size()==orb_A2.size());
	assert(orb_B1.size()==orb_B2.size());
	for( size_t TA = 0; TA!=orb_A1.size(); ++TA )
		for( size_t TB=0; TB!=orb_B1.size(); ++TB )
			for( int LA1=0; LA1!=orb_A1[TA].size(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1[TA][LA1].size(); ++NA1 )
					for( int LA2=0; LA2!=orb_A2[TA].size(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2[TA][LA2].size(); ++NA2 )
							for( int LB1=0; LB1!=orb_B1[TB].size(); ++LB1 )
								for( size_t NB1=0; NB1!=orb_B1[TB][LB1].size(); ++NB1 )
									for( size_t LB2=0; LB2!=orb_B2[TB].size(); ++LB2 )
										for( size_t NB2=0; NB2!=orb_B2[TB][LB2].size(); ++NB2 )
											center2_orb22_s[TA][TB][LA1][NA1][LA2][NA2][LB1][NB1][LB2].insert(
												make_pair(NB2, Center2_Orb::Orb22(
													orb_A1[TA][LA1][NA1],
													orb_A2[TA][LA2][NA2],
													orb_B1[TB][LB1][NB1],
													orb_B2[TB][LB2][NB2],
													MOT, MGT)));
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::init_radial\t"<<time_during(t_start)<<endl;
ofs.close();
}

void Exx_Abfs::Matrix_Orbs22::init_radial(
	const LCAO_Orbitals &orb_A1,
	const LCAO_Orbitals &orb_A2,
	const LCAO_Orbitals &orb_B1,
	const LCAO_Orbitals &orb_B2 )
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	TITLE("Exx_Abfs::Matrix_Orbs22","init_radial");

	assert( orb_A1.get_ntype() == orb_A2.get_ntype() );
	assert( orb_B1.get_ntype() == orb_B2.get_ntype() );
	for( size_t TA=0; TA!=orb_A1.get_ntype(); ++TA )
		for( size_t TB=0; TB!=orb_B1.get_ntype(); ++TB )
			for( int LA1=0; LA1<=orb_A1.Phi[TA].getLmax(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1.Phi[TA].getNchi(LA1); ++NA1 )
					for( int LA2=0; LA2<=orb_A2.Phi[TA].getLmax(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2.Phi[TA].getNchi(LA2); ++NA2 )
							for( int LB1=0; LB1<=orb_B1.Phi[TB].getLmax(); ++LB1 )
								for( size_t NB1=0; NB1!=orb_B1.Phi[TB].getNchi(LB1); ++NB1 )
									for( int LB2=0; LB2<=orb_B2.Phi[TB].getLmax(); ++LB2 )
										for( size_t NB2=0; NB2!=orb_B2.Phi[TB].getNchi(LB2); ++NB2 )
											center2_orb22_s[TA][TB][LA1][NA1][LA2][NA2][LB1][NB1][LB2].insert(
												make_pair(NB2, Center2_Orb::Orb22(
													orb_A1.Phi[TA].PhiLN(LA1,NA1),
													orb_A2.Phi[TA].PhiLN(LA2,NA2),
													orb_B1.Phi[TB].PhiLN(LB1,NB1),
													orb_B2.Phi[TB].PhiLN(LB2,NB2),
													MOT, MGT)));
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::init_radial\t"<<time_during(t_start)<<endl;
ofs.close();
}

void Exx_Abfs::Matrix_Orbs22::init_radial_table()
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	TITLE("Exx_Abfs::Matrix_Orbs22","init_radial_table");

	for( auto &coA : center2_orb22_s )
		for( auto &coB : coA.second )
			for( auto &coC : coB.second )
				for( auto &coD : coC.second )
					for( auto &coE : coD.second )
						for( auto &coF : coE.second )
							for( auto &coG : coF.second )
								for( auto &coH : coG.second )
									for( auto &coI : coH.second )
										for( auto &coJ : coI.second )
											coJ.second.init_radial_table();
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::init_radial_table\t"<<time_during(t_start)<<endl;
ofs.close();
}

void Exx_Abfs::Matrix_Orbs22::init_radial_table( const map<size_t,map<size_t,set<double>>> &Rs )
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	TITLE("Exx_Abfs::Matrix_Orbs22","init_radial_table_Rs");

	for( const auto &RsA : Rs )
		for( const auto &RsB : RsA.second )
		{
			if( auto* center2_orb22_sAB = static_cast<map<int,map<size_t,map<int,map<size_t,map<int,map<size_t,map<int,map<size_t,Center2_Orb::Orb22>>>>>>>>*>(
						MAP_EXIST(center2_orb22_s, RsA.first, RsB.first)) )
			{
timeval t_small;
gettimeofday(&t_small, NULL);
				set<size_t> radials;
				for( const double &R : RsB.second )
				{
					const double position = R * ucell.lat0 / MOT.dr;
					const size_t iq = static_cast<size_t>(position);
					for( size_t i=0; i!=4; ++i )
						radials.insert(iq+i);
				}
ofs<<"\t"<<RsA.first<<"\t"<<RsB.first<<"\t"<<time_during(t_small)<<"\t"<<flush;
gettimeofday(&t_small, NULL);
				for( auto &coC : *center2_orb22_sAB )
					for( auto &coD : coC.second )
						for( auto &coE : coD.second )
							for( auto &coF : coE.second )
								for( auto &coG : coF.second )
									for( auto &coH : coG.second )
										for( auto &coI : coH.second )
											for( auto &coJ : coI.second )
												coJ.second.init_radial_table(radials);
ofs<<time_during(t_small)<<endl;
			}
		}
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::init_radial_table_Rs\t"<<time_during(t_start)<<endl;
ofs.close();
}


matrix Exx_Abfs::Matrix_Orbs22::cal_overlap_matrix(
	const size_t TA,
	const size_t TB,
	const Vector3<double> &tauA,
	const Vector3<double> &tauB,
	const Element_Basis_Index::IndexLNM &index_A1,
	const Element_Basis_Index::IndexLNM &index_A2,
	const Element_Basis_Index::IndexLNM &index_B1,
	const Element_Basis_Index::IndexLNM &index_B2,
	const Matrix_Order &matrix_order) const
{
	TITLE("Exx_Abfs::Matrix_Orbs22","cal_overlap_matrix");

	matrix m;
	switch(matrix_order)
	{
		case Matrix_Order::A1B1_A2B2:	m.create( index_A1[TA].count_size*index_B1[TB].count_size, index_A2[TA].count_size*index_B2[TB].count_size );	break;
		default:	throw invalid_argument( "Matrix_Order wrong in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__) );
	}

	for( const auto &co3 : center2_orb22_s.at(TA).at(TB) )
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
								const int LB1 = co7.first;
								for( const auto &co8 : co7.second )
								{
									const size_t NB1 = co8.first;
									for( size_t MB1=0; MB1!=2*LB1+1; ++MB1 )
									{
										for( const auto &co9 : co8.second )
										{
											const int LB2 = co9.first;
											for( const auto & co10 : co9.second )
											{
												const size_t NB2 = co10.first;
												for( size_t MB2=0; MB2!=2*LB2+1; ++MB2 )
												{

													const double overlap = co10.second.cal_overlap( tauA*ucell.lat0, tauB*ucell.lat0, MA1, MA2, MB1, MB2 );

													switch(matrix_order)
													{
														case Matrix_Order::A1B1_A2B2:
															m(
																Exx_Abfs::Abfs_Index::get_index_index(
																	index_A1,TA,LA1,NA1,MA1,
																	index_B1,TB,LB1,NB1,MB1 ),
																Exx_Abfs::Abfs_Index::get_index_index(
																	index_A2,TA,LA2,NA2,MA2,
																	index_B2,TB,LB2,NB2,MB2 ) )
															= overlap;
															break;
														default:	throw invalid_argument( "Matrix_Order wrong in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__) );
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
	return m;
}


map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> Exx_Abfs::Matrix_Orbs22::cal_overlap_matrix(
	const Element_Basis_Index::IndexLNM &index_A1,
	const Element_Basis_Index::IndexLNM &index_A2,
	const Element_Basis_Index::IndexLNM &index_B1,
	const Element_Basis_Index::IndexLNM &index_B2 ) const
{
ofstream ofs(exx_lcao.test_dir.process+"time_"+TO_STRING(MY_RANK),ofstream::app);
timeval t_start;
gettimeofday( &t_start, NULL);
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> matrixes;

	for( const auto &co1 : center2_orb22_s )
	{
		const size_t TA = co1.first;
		for( size_t IA=0; IA!=ucell.atoms[TA].na; ++IA )
		{
			const Vector3<double> &tauA( ucell.atoms[TA].tau[IA] );

			for( const auto &co2 : co1.second )
			{
				const size_t TB = co2.first;
				for( size_t IB=0; IB!=ucell.atoms[TB].na; ++IB )
				{
					const Vector3<double> &tauB( ucell.atoms[TB].tau[IB] );

					matrixes[TA][IA][TB][IB] = cal_overlap_matrix(
						TA,
						TB,
						ucell.atoms[TA].tau[IA],
						ucell.atoms[TB].tau[IB],
						index_A1,
						index_A2,
						index_B1,
						index_B2,
						Matrix_Order::A1B1_A2B2);
				}
			}
		}
	}
ofs<<"TIME@Exx_Abfs::Matrix_Orbs22::cal_overlap_matrix\t"<<time_during(t_start)<<endl;
ofs.close();
	return matrixes;
}
