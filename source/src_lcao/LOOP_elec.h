#ifndef LOOP_ELEC_H
#define LOOP_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"

class LOOP_elec
{
	public:

	LOOP_elec(){};
	~LOOP_elec(){};

	// mohan add 2021-02-09
    void solve_elec_stru(const int& istep,
        Record_adj &ra,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& low,
        LCAO_Hamilt& uhm_in);

	private:

	// set matrix and grid integral
	void set_matrix_grid(Record_adj &ra);

    void before_solver(const int& istep,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf);

    void solver(const int& istep,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf);

    LCAO_Hamilt* UHM;

};

#endif
