#ifndef CENTER2_ORB_ORB22_CCP_H
#define CENTER2_ORB_ORB22_CCP_H

#include "center2_orb-orb22.h"

class Center2_Orb::Orb22_Ccp: public Center2_Orb::Orb22
{
public:
//	using Center2_Orb::Orb22::Orb22;
	Orb22_Ccp(
		const Numerical_Orbital_Lm &nA1_in,
		const Numerical_Orbital_Lm &nA2_in,
		const Numerical_Orbital_Lm &nB1_in,
		const Numerical_Orbital_Lm &nB2_in,		
		const ORB_table_phi &MOT_in,
		const ORB_gaunt_table &MGT_in	)
		: Orb22( nA1_in, nA2_in, nB1_in, nB2_in, MOT_in, MGT_in ){}

	void init_radial_table();
	Orb22_Ccp &set_rmesh_times( double rmesh_times_in ){ rmesh_times = rmesh_times_in; return *this; }
private:
	double rmesh_times = 1.0;
};

#endif