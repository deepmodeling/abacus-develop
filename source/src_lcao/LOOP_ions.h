#ifndef LOOP_IONS_H
#define LOOP_IONS_H

#include "LOOP_elec.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_pw/charge_extra.h"
#include "../src_ions/lattice_change_methods.h"
#include "src_lcao/local_orbital_wfc.h"
#include "module_orbital/ORB_control.h"
#include "src_lcao/LCAO_hamilt.h"
#include "module_esolver/esolver.h"

#include <fstream>

class LOOP_ions
{

	public:

	LOOP_ions(LCAO_Matrix &lm);
	~LOOP_ions();

    LOOP_elec LOE;
    Local_Orbital_wfc LOWF;
    Local_Orbital_Charge LOC;
    LCAO_Hamilt UHM;

	void opt_ions(ModuleESolver::ESolver *p_esolver); //output for dos
	void output_HS_R(
        const std::string &SR_filename="data-SR-sparse_SPIN0.csr",
        const std::string &HR_filename_up="data-HR-sparse_SPIN0.csr",
        const std::string HR_filename_down="data-HR-sparse_SPIN1.csr",
        const bool &binary=false, 
        const double &sparse_threshold=1e-10
    ); //LiuXh add 2019-07-15, modify in 2021-12-3
    void output_SR(const std::string &SR_filename, const bool &binary=false, const double &sparse_threshold=1e-10);

	private:

	Ions_Move_Methods IMM;

    Lattice_Change_Methods LCM;
	
	// PLEASE move 'force_stress()'  function to other places, such as FORCE_STRESS.cpp or
	// you might think to create a new file, it is because 'force_stress' do not
	// belong to 'LOOP_ions', 'GlobalC::pot.init_pot' also do not belong to force_stress()
	// the renew of structure factors, etc. should be ran in other places
	// the 'IMM' and 'LCM' objects should be passed to force_stress() via parameters list
	// mohan note 2021-03-23
	bool force_stress(const int &istep, int &force_step, int &stress_step, Record_adj &ra);

	int istep;

	// electron charge density extropolation method
	Charge_Extra CE;

	// PLEASE move final_scf to other places, for example, now it can be put 
	// in 'src_ions/variable_cell.cpp'
	// mohan note 2021-03-23
	void final_scf(void);
};

#endif
