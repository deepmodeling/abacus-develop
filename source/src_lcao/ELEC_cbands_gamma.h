#ifndef ELEC_CBANDS_GAMMA_H
#define ELEC_CBANDS_GAMMA_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "LCAO_hamilt.h"
#include "src_lcao/local_orbital_wfc.h"

//-----------------------------------------------------------
// mohan add 2021-02-09
// This class is used to calculate the band structures
// (eigenvalues and eigen wave functions) of the 
// Kohn-Sham equation in terms of the gamma-only k-point
// k is the index for the points in the first Brillouin zone
//-----------------------------------------------------------

class ELEC_cbands_gamma
{

	friend class ELEC_scf; 
	friend class ELEC_nscf; 

	public:

	ELEC_cbands_gamma();
	~ELEC_cbands_gamma();


	private:

    static void cal_bands(const int& istep, LCAO_Hamilt& uhm,
        Local_Orbital_wfc &lowf,
        std::vector<ModuleBase::matrix>& dm_gamma);


};

#endif
