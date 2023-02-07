#ifndef CONV_COULOMB_POT_UNITTEST_H
#define CONV_COULOMB_POT_UNITTEST_H

#include "src_ri/conv_coulomb_pot.h"
#include "test_function.h"

static void test_ccp()
{
	Conv_Coulomb_Pot ccp(GlobalC::ORB.Phi[0].PhiLN(0,0));
	ccp.cal_conv_coulomb_pot();
	std::stringstream ss;
	for( int ir=0; ir<=4*GlobalC::ORB.Phi[0].PhiLN(0,0).getNr(); ++ir )
	{
		ss<<ccp.get_conv_coulomb_pot(ir)<<std::endl;
	}
	MPI_RANK_OFSTREAM( "Conv_Coulomb_Pot", ss );
}

#endif