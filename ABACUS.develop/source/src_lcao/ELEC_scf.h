#ifndef ELEC_SCF_H
#define ELEC_SCF_H

#include "../src_pw/tools.h"
#include "../src_pw/threshold_elec.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to run self-consistent calculations 
// to solve the Kohn-Sham equation 
//-----------------------------------------------------------

class ELEC_scf: private Threshold_Elec
{

	friend class LOOP_elec;
	friend class LOOP_ions;
	friend class Run_MD;


	public:

	ELEC_scf();
	~ELEC_scf();

	private:

	void scf(const int &istep);

	static int iter;

	complex<double>*** WFC_init;

	void init_mixstep_final_scf(void);

};

#endif
