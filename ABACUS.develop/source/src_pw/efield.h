#ifndef EFIELD_H
#define EFIELD_H

#include "../module_base/matrix.h"

class Efield
{
	public:

	Efield();
	~Efield();

	void add_efield(const double*const rho, double* v_in);
    
	// efield.
    static int edir;         // direction of the field
    static double emaxpos;   // position of the maximum of the field (0<emaxpos<1)
    static double eopreg;    // amplitude of the inverse region (0<eopreg<1)
    static double eamp;    // field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
    static double etotefield;  // energy correction due to the field
	static double bvec[3];
	static double bmod;

	static void compute_force(matrix &fdip);

	private:

	void compute_el_dip(const double& emaxpos, const double &eopreg, 
		const int &edir, const double*const rho, double& e_dipole)const;

	void compute_ion_dip(const double& emaxpos, const double& eoperg, const int &edir, 
		double& ion_dipole, const double& bmod, const double* bg)const;

	double saw(const double &emaxpos, const double &eopreg, const double &x)const;


};

#endif
