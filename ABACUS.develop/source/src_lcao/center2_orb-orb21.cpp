//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-01-24
//=========================================================

#include "center2_orb-orb21.h"

#include "src_global/constants.h"
#include "src_global/ylm.h"

#include <cmath>

Center2_Orb::Orb21::Orb21(
	const Numerical_Orbital_Lm &nA1_in,
	const Numerical_Orbital_Lm &nA2_in,
	const Numerical_Orbital_Lm &nB_in,
	const ORB_table_phi &MOT_in,
	const ORB_gaunt_table &MGT_in
)
	:nA1(nA1_in),
	 nA2(nA2_in),
	 nB(nB_in),
	 MOT(MOT_in),
	 MGT(MGT_in)
{}

void Center2_Orb::Orb21::init_radial_table()
{
	const Numerical_Orbital_Lm & nA_short = (nA1.getNr()<=nA2.getNr()) ? nA1 : nA2;
	
	vector<double> nA_tmp( nA_short.getNr() );
	for( size_t ir=0; ir!=nA_tmp.size(); ++ir)
	{
		nA_tmp[ir] = nA1.getPsi(ir) * nA2.getPsi(ir);
	}
	
	const int LA1 = nA1.getL();
	const int LA2 = nA2.getL();
	for( int LA = std::abs(LA1-LA2); LA<=LA1+LA2; ++LA)
	{
		if( (LA-std::abs(LA1-LA2))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;		

		nA[LA].set_orbital_info(
			nA_short.getLabel(),
			nA_short.getType(),
			LA,
			1,						// N?
			nA_short.getNr(),
			nA_short.getRab(),
			nA_short.getRadial(),
			Numerical_Orbital_Lm::Psi_Type::Psi,
			VECTOR_TO_PTR(nA_tmp),
			nA_short.getNk(),
			nA_short.getDk(),
			nA_short.getDruniform(),
			false,
			true,
			FORCE); // mohan add 2021-05-07

		orb11s.insert( make_pair( LA, Center2_Orb::Orb11(nA[LA], nB, MOT, MGT) ) );

		orb11s.at(LA).init_radial_table();

	}
}

void Center2_Orb::Orb21::init_radial_table( const set<size_t> &radials )
{
	const Numerical_Orbital_Lm & nA_short = (nA1.getNr()<=nA2.getNr()) ? nA1 : nA2;
	
	vector<double> nA_tmp( nA_short.getNr() );
	for( size_t ir=0; ir!=nA_tmp.size(); ++ir)
	{
		nA_tmp[ir] = nA1.getPsi(ir) * nA2.getPsi(ir);
	}
	
	const int LA1 = nA1.getL();
	const int LA2 = nA2.getL();
	for( int LA = std::abs(LA1-LA2); LA<=LA1+LA2; ++LA)
	{
		if( (LA-std::abs(LA1-LA2))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;		

		nA[LA].set_orbital_info(
			nA_short.getLabel(),
			nA_short.getType(),
			LA,
			1,						// N?
			nA_short.getNr(),
			nA_short.getRab(),
			nA_short.getRadial(),
			Numerical_Orbital_Lm::Psi_Type::Psi,
			VECTOR_TO_PTR(nA_tmp),
			nA_short.getNk(),
			nA_short.getDk(),
			nA_short.getDruniform(),
			false,
			true, 
			FORCE); // mohan add 2021-05-07

		orb11s.insert( make_pair( LA, Center2_Orb::Orb11(nA[LA], nB, MOT, MGT) ) );

		orb11s.at(LA).init_radial_table(radials);

	}
}

double Center2_Orb::Orb21::cal_overlap(
	const Vector3<double> &RA, const Vector3<double> &RB,
	const int &mA1, const int &mA2, const int &mB) const
{
	const int LA1 = nA1.getL();
	const int LA2 = nA2.getL();
	
	double overlap = 0.0;
	
	for( const auto& orb11 : orb11s )
	{
		const int LA = orb11.first;

		for( int mA=0; mA!=2*LA+1; ++mA )
		// const int mA=mA1+mA2;
		{
			const double Gaunt_real_A1_A2_A12 = 
				MGT.Gaunt_Coefficients (
					MGT.get_lm_index(LA1,mA1),
					MGT.get_lm_index(LA2,mA2),
					MGT.get_lm_index(LA,mA));
			if( 0==Gaunt_real_A1_A2_A12 )	continue;

			overlap += Gaunt_real_A1_A2_A12 * orb11.second.cal_overlap(RA, RB, mA, mB);
		}
	}
	
	return overlap;
}
