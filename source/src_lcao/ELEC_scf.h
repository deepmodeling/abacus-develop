#ifndef ELEC_SCF_H
#define ELEC_SCF_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "../src_pw/threshold_elec.h"
#include "src_lcao/local_orbital_charge.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to run self-consistent calculations 
// to solve the Kohn-Sham equation 
//-----------------------------------------------------------

class ELEC_scf: private Threshold_Elec
{

	friend class LOOP_elec;
	friend class LOOP_ions;
	friend class Run_MD_LCAO;


	public:

	ELEC_scf();
	~ELEC_scf();

	//private:

    void scf(const int& istep,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf,
        LCAO_Hamilt& uhm);

	static int iter;

	void init_mixstep_final_scf(void);

};

#endif
