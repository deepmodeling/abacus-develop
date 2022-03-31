#ifndef ELEC_NSCF_H
#define ELEC_NSCF_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "LCAO_hamilt.h"
#include "src_lcao/local_orbital_wfc.h"
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

    static void nscf(LCAO_Hamilt& uhm,
        std::vector<ModuleBase::matrix>& dm_gamma,
        std::vector<ModuleBase::ComplexMatrix>& dm_k,
        Local_Orbital_wfc &lowf);


};

#endif
