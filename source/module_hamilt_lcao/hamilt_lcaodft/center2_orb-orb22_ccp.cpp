#include "center2_orb-orb22_ccp.h"
#include "conv_coulomb_pot.h"
#include "conv_coulomb_pot-inl.h"

void Center2_Orb::Orb22_Ccp::init_radial_table()
{
	std::vector<double> nB_tmp(nB1.getNr());
	for( size_t ir=0; ir!=nB_tmp.size(); ++ir)
	{
		nB_tmp[ir] = nB1.getPsi(ir) * nB2.getPsi(ir);
	}	
	const int LB1 = nB1.getL();
	const int LB2 = nB2.getL();
	for( int LB = abs(LB1-LB2); LB<=LB1+LB2; ++LB)
	{
		if( (LB-std::abs(LB1-LB2))%2==1 )			// if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
			continue;

//std::cout<<"LB\t"<<LB<<std::endl;
			
		Numerical_Orbital_Lm nB_origin;
		nB_origin.set_orbital_info(
			nB1.getLabel(),
			nB1.getType(),
			LB,
			1,						// N?
			nB1.getNr(),
			nB1.getRab(),
			nB1.getRadial(),
			ModuleBase::GlobalFunc::VECTOR_TO_PTR(nB_tmp),
			nB1.getNk(),
			nB1.getDk(),
			nB1.getDruniform(),
			false,
			true);

//std::cout<<"before ccp"<<std::endl;
		Conv_Coulomb_Pot::cal_orbs_ccp( nB_origin, nB[LB], rmesh_times, 1 );		// Peize Lin test
//		nB[LB].cal_kradial();
//std::cout<<"after ccp"<<std::endl;

		orb21s.insert( make_pair( LB, Center2_Orb::Orb21( nA1, nA2, nB[LB], MOT, MGT ) ) );

//std::cout<<"before orb21s"<<std::endl;		
		orb21s.at(LB).init_radial_table();
//std::cout<<"after orb21s"<<std::endl;		
	}
}