#ifndef NUMERICAL_ORBITAL_H
#define NUMERICAL_ORBITAL_H

#include <string>
using namespace std;

#include "../module_base/intarray.h"
#include "../module_base/vector3.h"
#include "ORB_atomic_lm.h"

///
///CLASS  Num_Orbital
///------------------------------------------
///
///Note : contain information about atoms
///
///Feature : set and store information about atoms
///
class Numerical_Orbital
{
	friend class LCAO_Orbitals;

public:

	Numerical_Orbital();
	~Numerical_Orbital();

	const int& getLmax() const { return this->lmax; }
	const double& getRcut () const {return this->rcut; }
   	const int& getType() const { return this->type; }
	const int& getTotal_nchi() const { return this->total_nchi; }
	const int& getNchi(const int l) const { return this->nchi[l]; }
	const string& getLabel() const { return this->label; }
	
	const inline Numerical_Orbital_Lm& PhiLN( const int &L, const int &N)const
	{ 	
		return this->phiLN[ this->find_chi(L, N) ];
	}
	
	/// about the distance between two atoms.
	static double& get_distance()
	{ 
		if(distance < 0.0) WARNING_QUIT("NUMERICAL_ORBITAL","distance should be above zero!"); 
		return distance; 
	}
	
	static double getX() { return R2.x - R1.x ; }
	static double getY() { return R2.y - R1.y ; }
	static double getZ() { return R2.z - R1.z ; }
	static Vector3<double>& getR1() { return R1; }
	static Vector3<double>& getR2() { return R2; }
	static Vector3<double>& getdR() { return dR; }

	///
	/// set information about Numerical Orbital
	///
	void set_orbital_info(
			const int& type_in,
			const string& label_in,
			const int& lmax_in,
			const int* nchi_in,
			const int& total_nchi);

	static void set_position(const Vector3<double>R1_in, const Vector3<double> R2_in) 
	{
		R1 = R1_in;
		R2 = R2_in;
		dR = R1-R2;
		distance = dR.norm();
	}
				
private:
	
	//==========================================================
	// MEMBER FUNCTION :
	// NAME : label (atom type)
	// NAME : lmax( max value of L angular momentum) 
	// NAME : nchi( number of chi for each L)
	// NAME : total_nchi(total chi for this type of atom, total number of NAOs)
	// NAME : max_nchi( max chi for certain L)
	// NAME : find_chi(lmax+1, max_nchi).
	//===========================================================
	string label;
	
	int type;
	int lmax;
	int* nchi;
	int total_nchi;
	int max_nchi;
	IntArray find_chi;
	double rcut;

	Numerical_Orbital_Lm* phiLN;// length: total_nchi (only store radial function )

	//==========================================================
	// It's about two atoms, so here we set static variables 
	//==========================================================
	static double distance; 
	static Vector3<double> R1;
	static Vector3<double> R2; //three-dimesion-coordinate of R
	static Vector3<double> dR; // R1-R2
};

#endif
