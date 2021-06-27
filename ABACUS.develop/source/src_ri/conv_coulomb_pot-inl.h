#ifndef CONV_COULOMB_POT_INL_H
#define CONV_COULOMB_POT_INL_H

#include "conv_coulomb_pot.h"
#include "../module_base/global_function.h"

#include "../src_external/src_test/test_function.h"			// Peize Lin test

inline double Conv_Coulomb_Pot::get_r_radial_outofrange( size_t ir ) const
{
	return orb.getRadial(conv_coulomb_pot.size()-1) 
		+ orb.getRab(conv_coulomb_pot.size()-1) 
		* ( ir - conv_coulomb_pot.size() + 1 ) ;	
}
	
template< typename T >
void Conv_Coulomb_Pot::cal_orbs_ccp( 
	const T & orbs, 
	T & orbs_ccp, 
	size_t rmesh_times, 
	size_t kmesh_times )
{
	orbs_ccp.resize(orbs.size());
	for( int i=0; i!=orbs.size(); ++i )
	{
		cal_orbs_ccp( orbs[i], orbs_ccp[i], rmesh_times, kmesh_times );
	}
}

template<>
void Conv_Coulomb_Pot::cal_orbs_ccp<Numerical_Orbital_Lm>( 
	const Numerical_Orbital_Lm & orbs, 
	Numerical_Orbital_Lm & orbs_ccp, 
	size_t rmesh_times, 
	size_t kmesh_times );

#endif
