#ifndef ELEC_NSCF_H
#define ELEC_NSCF_H

#include "../src_pw/tools.h"
#include "LCAO_hamilt.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to run non-self-consistent calculations 
// to solve the Kohn-Sham equation by reading in electron
// charge density
//-----------------------------------------------------------

class ELEC_nscf
{

	friend class LOOP_ions;
	friend class LOOP_elec;

	public:

	ELEC_nscf();
	~ELEC_nscf();

	private:

	static void nscf(LCAO_Hamilt &uhm); 


};

#endif
