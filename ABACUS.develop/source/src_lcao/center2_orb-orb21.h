//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-01-24
//=========================================================

#ifndef CENTER2_ORB_ORB21_H
#define CENTER2_ORB_ORB21_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "center2_orb.h"
#include "center2_orb-orb11.h"

#include "src_lcao/ORB_table_phi.h"
#include "src_lcao/ORB_gaunt_table.h"
#include "src_lcao/ORB_atomic_lm.h"
#include "src_global/vector3.h"

class Center2_Orb::Orb21
{

	public:

	Orb21(
		const Numerical_Orbital_Lm &nA1_in,
		const Numerical_Orbital_Lm &nA2_in,
		const Numerical_Orbital_Lm &nB_in,
		const ORB_table_phi &MOT_in,
		const ORB_gaunt_table &MGT_in	);

	void init_radial_table();

	void init_radial_table( const set<size_t> &radials );					// unit: Bohr/MOT.dr

	double cal_overlap(
		const Vector3<double> &RA, const Vector3<double> &RB,				// unit: Bohr
		const int &mA1, const int &mA2, const int &mB ) const;
		
	private:

	const Numerical_Orbital_Lm &nA1;
	const Numerical_Orbital_Lm &nA2;
	const Numerical_Orbital_Lm &nB;
	
	const ORB_table_phi &MOT;
	const ORB_gaunt_table &MGT;	

	map<int,Numerical_Orbital_Lm> nA;
	map<int,Center2_Orb::Orb11> orb11s;
};

// this->orb11s[LA].Table_r[LAB][ir]

#endif	// CENTER2_ORB_ORB21_H
