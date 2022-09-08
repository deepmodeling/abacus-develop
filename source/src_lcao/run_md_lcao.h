#ifndef RUN_MD_LCAO_H
#define RUN_MD_LCAO_H 

#include "../src_pw/charge_extra.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_matrix.h"
#include "module_esolver/esolver.h"

class Run_MD_LCAO
{

	public:

	Run_MD_LCAO();
	~Run_MD_LCAO();

	void opt_ions(ModuleESolver::ESolver *p_esolver);
	void md_force_virial(ModuleESolver::ESolver *p_esolver,
		const int &istep,
        const int& numIon, 
        double &potential, 
        ModuleBase::Vector3<double>* force, 
        ModuleBase::matrix& virial);

	private:

	// electron charge density extropolation method
	Charge_Extra CE;
    bool cellchange;
};

#endif
