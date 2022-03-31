#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "MD_parameters.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../module_base/matrix.h"
#ifdef __LCAO
#include "../module_orbital/parallel_orbitals.h"
#endif

#include "module_esolver/esolver.h"
class MD_func
{
    public:

    MD_func(){};
    ~MD_func() {};
    
#ifdef __LCAO
    static Parallel_Orbitals* ParaV;
#endif

    static double gaussrand();
	
	static void InitPos( const UnitCell_pseudo &unit_in, ModuleBase::Vector3<double>* pos);

	static void InitVel(
		const UnitCell_pseudo &unit_in, 
		const double& temperature, 
		double* allmass,
		int& frozen_freedom,
		ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	static void ReadVel(const UnitCell_pseudo &unit_in, ModuleBase::Vector3<double>* vel);

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
		const MD_parameters &mdp,
		const UnitCell_pseudo &unit_in,
		double &potential,
		ModuleBase::Vector3<double> *force,
		ModuleBase::matrix &stress);

	static double GetAtomKE(
		const int &numIon,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass);

	static void kinetic_stress(
		const UnitCell_pseudo &unit_in,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass, 
		double &kinetic,
		ModuleBase::matrix &stress);

	static void outStress(const ModuleBase::matrix &virial, const ModuleBase::matrix &stress);

	static void MDdump(
		const int &step, 
		const UnitCell_pseudo &unit_in,
		const ModuleBase::matrix &virial, 
		const ModuleBase::Vector3<double> *force);

	static void getMassMbl(const UnitCell_pseudo &unit_in, 
		double* allmass, 
		ModuleBase::Vector3<int> &frozen,
		ModuleBase::Vector3<int>* ionmbl);
};
#endif
