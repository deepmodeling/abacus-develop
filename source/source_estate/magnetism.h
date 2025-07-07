#ifndef MAGNETISM_H
#define MAGNETISM_H

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/vector3.h"

class Magnetism
{
public:
    // constructor and deconstructor
    Magnetism();
    ~Magnetism();

    // notice : bcast (MPI operation) is done in unitcell
    double *start_mag=nullptr;

    // tot_mag : majority spin - minority spin (nelup - neldw).
    double tot_mag;

    double tot_mag_nc[3]={0.0};

    double abs_mag;

	void compute_mag(const double& omega,
			const int& nrxx, 
			const int& nxyz, 
			const double* const * rho, 
			double* nelec_spin = nullptr);

    ModuleBase::Vector3<double> *m_loc_=nullptr; //magnetization for each element along c-axis

	double *angle1_=nullptr;                     //angle between c-axis and real spin std::vector

	double *angle2_=nullptr;                     //angle between a-axis and real spin std::vector projection in ab-plane

    double ux_[3]={0.0};

	bool lsign_=false;

private:

    bool judge_parallel(const double a[3], const ModuleBase::Vector3<double> &b);

};

/*
 A comment about variables nelup, neldw, multiplicity and tot_mag:
 All these variables contain the same information and must be kept harmonized.
 Variables nelup and neldw will be removed in future versions of the code.
 Variables multiplicity and tot_mag, though redundent will probably
 coexist since multiplicity is the more natural way (?)for defining the spin
 configuratio in the quantum-chemistry community while tot_mag is
 more natural (?) when dealing with extended systems.
*/

#endif
