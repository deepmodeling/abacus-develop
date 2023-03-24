#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "MD_parameters.h"
#include "../module_cell/unitcell.h"
#include "../module_base/matrix.h"

#include "module_esolver/esolver.h"
class MD_func
{
    public:

    MD_func(){};
    ~MD_func() {};

    static double gaussrand();
	
	static void InitPos( const UnitCell &unit_in, ModuleBase::Vector3<double>* pos);

	static void InitVel(
		const UnitCell &unit_in, 
		const double& temperature, 
		double* allmass,
		int& frozen_freedom,
		ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	static void ReadVel(const UnitCell &unit_in, ModuleBase::Vector3<double>* vel);

	static void RandomVel(
		const int& numIon, 
		const double& temperature,
		const double* allmass,
		const int& frozen_freedom,
		const ModuleBase::Vector3<int> frozen,
		const ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	static void force_virial(
		ModuleESolver::ESolver *p_esolver,
		const int &istep,
		UnitCell &unit_in,
		double &potential,
		ModuleBase::Vector3<double> *force,
		ModuleBase::matrix &virial);

	static double GetAtomKE(
		const int &numIon,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass);

	static void compute_stress(
		const UnitCell &unit_in,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass, 
        const ModuleBase::matrix &virial,
		ModuleBase::matrix &stress);

	static void outStress(const ModuleBase::matrix &virial, const ModuleBase::matrix &stress);

    static void MDdump(
        const int &step, 
        const UnitCell &unit_in,
        const Input &inp,
        const ModuleBase::matrix &virial, 
        const ModuleBase::Vector3<double> *force,
        const ModuleBase::Vector3<double> *vel);

	static void getMassMbl(const UnitCell &unit_in, 
		double* allmass, 
		ModuleBase::Vector3<int> &frozen,
		ModuleBase::Vector3<int>* ionmbl);

    static double target_temp(const int &istep, const double &tfirst, const double &tlast);

    static double current_temp(double &kinetic,
            const int &natom, 
            const int &frozen_freedom, 
            const double *allmass,
            const ModuleBase::Vector3<double> *vel);

    static void temp_vector(const int &natom, 
            const ModuleBase::Vector3<double> *vel, 
            const double *allmass, 
            ModuleBase::matrix &t_vector);
};
#endif
