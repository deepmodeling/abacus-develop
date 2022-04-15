#ifndef ELEC_EVOLVE_H
#define ELEC_EVOLVE_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "LCAO_hamilt.h"
#include "src_lcao/local_orbital_wfc.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to evolve the electronic wave functions
// in TDDFT in terms of the multiple k points
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class ELEC_evolve
{

	friend class ELEC_scf;

	public:

	ELEC_evolve();
	~ELEC_evolve();

	// fuxiang add 2021-05-25

    static int tddft;
    static double td_scf_thr;
    static double td_dt;
    static double td_force_dt;
    static int td_val_elec_01;
    static int td_val_elec_02;
    static int td_val_elec_03;
    static int td_vext;
    static int td_vext_dire;
	static double td_timescale;
	static int td_vexttype;
	static int td_vextout;
	static int td_dipoleout;


	private:

    static void evolve_psi(const int& istep, LCAO_Hamilt& uhm,
    Local_Orbital_wfc &lowf);

};

#endif
