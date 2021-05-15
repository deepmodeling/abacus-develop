#ifndef CONV_COULOMB_POT_H
#define CONV_COULOMB_POT_H

#include <vector>
using std::vector;

#include "../module_ORB/ORB_atomic_lm.h"

// out of conv_coulomb_pot.size(), imagine all rab(i.e. dr) are equal to rab[conv_coulomb_pot.size()-1]
class Conv_Coulomb_Pot
{
public:
	Conv_Coulomb_Pot(const Numerical_Orbital_Lm &orb_in);

	const vector<double> &get_conv_coulomb_pot() const { return conv_coulomb_pot; }
	double get_conv_coulomb_pot(const size_t &ir) const;

	void cal_conv_coulomb_pot();

public:
	// only for <vector ... <vector <Numerical_Orbital_Lm> > ... >
	// the number of "vector" can be 0 to any
	template< typename T >
	static void cal_orbs_ccp( const T & orb, T & orb_ccp, size_t rmesh_times=1, size_t kmesh_times=1 );

private:
	const Numerical_Orbital_Lm &orb;
	vector<double> conv_coulomb_pot;
	double conv_coulomb_pot_extra;

	inline double get_r_radial_outofrange( size_t ir ) const;
};

#endif
