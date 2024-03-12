//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-01-24
//=========================================================

#ifndef CENTER2_ORB_ORB11_H
#define CENTER2_ORB_ORB11_H

#include <vector>
#include <map>
#include <set>

#include "center2_orb.h"

#include "module_basis/module_ao/ORB_table_phi.h"
#include "module_basis/module_ao/ORB_gaunt_table.h"
#include "module_basis/module_ao/ORB_atomic_lm.h"
#include "module_base/vector3.h"

class Center2_Orb::Orb11
{
public:

	Orb11(
		const Numerical_Orbital_Lm &nA_in,
		const Numerical_Orbital_Lm &nB_in,
		const ORB_table_phi &MOT_in,
		const ORB_gaunt_table &MGT_in);

	void init_radial_table(void);
	void init_radial_table( const std::set<size_t> &radials );				// unit: Bohr/MOT.dr

	double cal_overlap(
		const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,			// unit: Bohr
		const int &mA, const int &mB) const;

    ModuleBase::Vector3<double> cal_grad_overlap( //caoyu add 2021-11-19
		const ModuleBase::Vector3<double> &RA, const ModuleBase::Vector3<double> &RB,			// unit: Bohr
		const int& mA, const int& mB) const;
    
private:

	const Numerical_Orbital_Lm &nA;
	const Numerical_Orbital_Lm &nB;
	
	const ORB_table_phi &MOT;
	const ORB_gaunt_table &MGT;	

	std::map<int,std::vector<double>> Table_r;									// unit: Bohr/MOT.dr
	std::map<int,std::vector<double>> Table_dr;
};

// this->Table_r[LAB][ir]

#endif	// CENTER2_ORB_ORB11_H
