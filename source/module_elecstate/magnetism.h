#ifndef MAGNETISM_H
#define MAGNETISM_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/vector3.h"
#include "module_elecstate/module_charge/charge.h"

class Magnetism
{
public:
    // constructor and deconstructor
    Magnetism();
    ~Magnetism();

    // notice : becast is done in unitcell
    double *start_magnetization;

    // tot_magnetization : majority spin - minority spin (nelup - neldw).
    double tot_magnetization;
    double tot_magnetization_nc[3];
    double abs_magnetization;

    void compute_magnetization(const int& nrxx, const int& nxyz, const double* const * rho, double* nelec_spin = nullptr);

    ModuleBase::Vector3<double> *m_loc_;   //magnetization for each element along c-axis
	double *angle1_;           //angle between c-axis and real spin std::vector
	double *angle2_;           //angle between a-axis and real spin std::vector projection in ab-plane
    //void cal_ux(const int ntype);
    double ux_[3];
	bool lsign_;

private:
    bool judge_parallel(double a[3],ModuleBase::Vector3<double> b);

};

/*
 A comment about variables nelup, neldw, multiplicity and tot_magnetization:
 All these variables contain the same information and must be kept harmonized.
 Variables nelup and neldw will be removed in future versions of the code.
 Variables multiplicity and tot_magnetization, though redundent will probably
 coexist since multiplicity is the more natural way (?)for defining the spin
 configuratio in the quantum-chemistry community while tot_magnetization is
 more natural (?) when dealing with extended systems.
*/

#endif
