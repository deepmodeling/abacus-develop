#ifndef NUMERICAL_ORBITAL_H
#define NUMERICAL_ORBITAL_H

#include <string>

#include "module_base/intarray.h"
#include "module_base/vector3.h"
#include "ORB_atomic_lm.h"

class Numerical_Orbital_AtomRelation
{
public:
	//==========================================================
	// It's about two atoms relations, thread-safe interface
	//==========================================================
	double distance; 
	ModuleBase::Vector3<double> R1;
	ModuleBase::Vector3<double> R2; //three-dimesion-coordinate of R
	ModuleBase::Vector3<double> dR; // R1-R2

	double& get_distance()
	{ 
		if(distance < 0.0) ModuleBase::WARNING_QUIT("NUMERICAL_ORBITAL","distance should be above zero!"); 
		return distance; 
	}
	
	double getX() { return R2.x - R1.x ; }
	double getY() { return R2.y - R1.y ; }
	double getZ() { return R2.z - R1.z ; }
	ModuleBase::Vector3<double>& getR1() { return R1; }
	ModuleBase::Vector3<double>& getR2() { return R2; }
	ModuleBase::Vector3<double>& getdR() { return dR; }

	void set_position(const ModuleBase::Vector3<double> &R1_in, const ModuleBase::Vector3<double> &R2_in) 
	{
		R1 = R1_in;
		R2 = R2_in;
		dR = R1-R2;
		distance = dR.norm();
	}
};

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
	const std::string& getLabel() const { return this->label; }
	
	const inline Numerical_Orbital_Lm& PhiLN( const int &L, const int &N)const
	{ 	
		return this->phiLN[ this->find_chi(L, N) ];
	}
	
	/// about the distance between two atoms.
	static double& get_distance()
	{
		return NOAR.get_distance(); 
	}
	
	static double getX() { return NOAR.getX() ; }
	static double getY() { return NOAR.getY() ; }
	static double getZ() { return NOAR.getZ() ; }
	static ModuleBase::Vector3<double>& getR1() { return NOAR.getR1(); }
	static ModuleBase::Vector3<double>& getR2() { return NOAR.getR2(); }
	static ModuleBase::Vector3<double>& getdR() { return NOAR.getdR(); }

	///
	/// set information about Numerical Orbital
	///
	void set_orbital_info(
			const int& type_in,
			const std::string& label_in,
			const int& lmax_in,
			const int* nchi_in,
			const int& total_nchi);

	static void set_position(const ModuleBase::Vector3<double> &R1_in, const ModuleBase::Vector3<double> &R2_in) 
	{
		NOAR.set_position(R1_in, R2_in);
	}

	Numerical_Orbital_Lm*& chi() { return this->phiLN; }
				
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
	std::string label;
	
	int type;
	int lmax;
	int* nchi;
	int total_nchi;
	int max_nchi;
	ModuleBase::IntArray find_chi;
	double rcut;

	Numerical_Orbital_Lm* phiLN;// length: total_nchi (only store radial function )

	//==========================================================
	// Keep the old interface
	//==========================================================
	static Numerical_Orbital_AtomRelation NOAR;
};

#endif
