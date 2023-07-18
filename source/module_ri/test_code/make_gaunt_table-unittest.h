//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-09-09
//==========================================================

#ifndef MAKE_GAUNT_TABLE_UNITTEST_H
#define MAKE_GAUNT_TABLE_UNITTEST_H

#include "module_basis/module_ao/ORB_gaunt_table.h"

static void cout_MGT ( const ORB_gaunt_table & MGT, const int Lmax )
{
	for( int LA=0; LA<=Lmax; ++LA )
	{
		for( int LB=0; LB<=Lmax; ++LB )
		{
			for( int LAB=std::abs(LA-LB); LAB<=LA+LB; ++LAB )
			{
				for( int mA=0; mA!=2*LA+1; ++mA )
				{
					for( int mB=0; mB!=2*LB+1; ++mB )
					{
						for( int mAB=0; mAB!=2*LAB+1; ++mAB )
						{
							std::cout<<LA<<"\t"
								<<LB<<"\t"
								<<LAB<<"\t"
								<<mA<<"\t"
								<<mB<<"\t"
								<<mAB<<"\t"
								<<ORB_gaunt_table::Index_M(mA)<<"\t"
								<<ORB_gaunt_table::Index_M(mB)<<"\t"
								<<ORB_gaunt_table::Index_M(mAB)<<"\t"
								<<MGT.Gaunt_Coefficients (
									MGT.get_lm_index(LA,mA),
									MGT.get_lm_index(LB,mB),
									MGT.get_lm_index(LAB,mAB))<<std::endl;
						}
					}
				}
			}

		}
	}
}
#endif
