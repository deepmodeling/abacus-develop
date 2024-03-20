//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#include "Matrix_Orbs21.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"

void Matrix_Orbs21::init(
	const int mode,
	const double kmesh_times,
	const double rmesh_times)
{
	ModuleBase::TITLE("Matrix_Orbs21","init");
	ModuleBase::timer::tick("Matrix_Orbs21", "init");
	//=========================================
	// (1) MOT: make overlap table.
	//=========================================
	this->MOT.allocate(
		GlobalC::ORB.get_ntype(),							// number of atom types
		std::max( GlobalC::ORB.get_lmax(), GlobalC::exx_info.info_ri.abfs_Lmax ),	// max L used to calculate overlap
		static_cast<int>(GlobalC::ORB.get_kmesh() * kmesh_times) | 1,				// kpoints, for integration in k space
		GlobalC::ORB.get_Rmax() * rmesh_times,				// max value of radial table
		GlobalC::ORB.get_dR(),								// delta R, for making radial table
//		GlobalC::ORB.get_dk() / kmesh_times);				// delta k, for integration in k space
		GlobalC::ORB.get_dk());											// Peize Lin change 2017-04-16
	int Lmax_used, Lmax;
	this->MOT.init_Table_Spherical_Bessel (3,mode, Lmax_used, Lmax, GlobalC::exx_info.info_ri.abfs_Lmax, GlobalC::ORB, GlobalC::ucell.infoNL.Beta);
//	this->MOT.init_OV_Tpair();							// for this->MOT.OV_L2plus1
//	this->MOT.Destroy_Table_Spherical_Bessel (Lmax_used);				// why?

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	ModuleBase::Ylm::set_coefficients ();

	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================
	this->MGT.init_Gaunt_CH( 2*Lmax+1 );			// why +1
	this->MGT.init_Gaunt( 2*Lmax+1 );

	ModuleBase::timer::tick("Matrix_Orbs21", "init");
}



void Matrix_Orbs21::init_radial(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A2,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_B )
{
	ModuleBase::TITLE("Matrix_Orbs21","init_radial");
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
	assert(orb_A1.size()==orb_A2.size());
	for( size_t TA=0; TA!=orb_A1.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.size(); ++TB )
			for( int LA1=0; LA1!=orb_A1[TA].size(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1[TA][LA1].size(); ++NA1 )
					for( int LA2=0; LA2!=orb_A2[TA].size(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2[TA][LA2].size(); ++NA2 )
							for( int LB=0; LB!=orb_B[TB].size(); ++LB )
								for( size_t NB=0; NB!=orb_B[TB][LB].size(); ++NB )
									center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
										std::make_pair(NB, Center2_Orb::Orb21(
											orb_A1[TA][LA1][NA1],
											orb_A2[TA][LA2][NA2],
											orb_B[TB][LB][NB],
											this->MOT, this->MGT)));
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
}


void Matrix_Orbs21::init_radial(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &orb_A1,
	const LCAO_Orbitals &orb_A2,
	const LCAO_Orbitals &orb_B)
{
	ModuleBase::TITLE("Matrix_Orbs21","init_radial");
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
	assert( orb_A1.size() == orb_A2.get_ntype() );
	for( size_t TA=0; TA!=orb_A1.size(); ++TA )
		for( size_t TB=0; TB!=orb_B.get_ntype(); ++TB)
			for( int LA1=0; LA1!=orb_A1[TA].size(); ++LA1 )
				for( size_t NA1=0; NA1!=orb_A1[TA][LA1].size(); ++NA1 )
					for( int LA2=0; LA2<=orb_A2.Phi[TA].getLmax(); ++LA2 )
						for( size_t NA2=0; NA2!=orb_A2.Phi[TA].getNchi(LA2); ++NA2 )
							for( int LB=0; LB<=orb_B.Phi[TB].getLmax(); ++LB )
								for( size_t NB=0; NB!=orb_B.Phi[TB].getNchi(LB); ++NB )
									center2_orb21_s[TA][TB][LA1][NA1][LA2][NA2][LB].insert(
										std::make_pair(NB, Center2_Orb::Orb21(
											orb_A1[TA][LA1][NA1],
											orb_A2.Phi[TA].PhiLN(LA2,NA2),
											orb_B.Phi[TB].PhiLN(LB,NB),
											this->MOT, this->MGT)));
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial");
}


void Matrix_Orbs21::init_radial_table()
{
	ModuleBase::TITLE("Matrix_Orbs21","init_radial");
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
	for( auto &coA : center2_orb21_s )
		for( auto &coB : coA.second )
			for( auto &coC : coB.second )
				for( auto &coD : coC.second )
					for( auto &coE : coD.second )
						for( auto &coF : coE.second )
							for( auto &coG : coF.second )
								for( auto &coH : coG.second )
									coH.second.init_radial_table();
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
}

void Matrix_Orbs21::init_radial_table( const std::map<size_t,std::map<size_t,std::set<double>>> &Rs )
{
	ModuleBase::TITLE("Matrix_Orbs21","init_radial_table_Rs");
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
	for( const auto &RsA : Rs )
		for( const auto &RsB : RsA.second )
		{
			if( auto* const center2_orb21_sAB = static_cast<std::map<int,std::map<size_t,std::map<int,std::map<size_t,std::map<int,std::map<size_t,Center2_Orb::Orb21>>>>>>*const>(
						ModuleBase::GlobalFunc::MAP_EXIST(center2_orb21_s, RsA.first, RsB.first)) )
			{
				std::set<size_t> radials;
				for( const double &R : RsB.second )
				{
					const double position = R * GlobalC::ucell.lat0 / this->MOT.dr;
					const size_t iq = static_cast<size_t>(position);
					for( size_t i=0; i!=4; ++i )
						radials.insert(iq+i);
				}
				for( auto &coC : *center2_orb21_sAB )
					for( auto &coD : coC.second )
						for( auto &coE : coD.second )
							for( auto &coF : coE.second )
								for( auto &coG : coF.second )
									for( auto &coH : coG.second )
										coH.second.init_radial_table(radials);
			}
		}
	ModuleBase::timer::tick("Matrix_Orbs21", "init_radial_table");
}
