#ifndef LCAO_SET_H
#define LCAO_SET_H

#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_psi/psi.h"
#include "source_estate/elecstate.h"
#include "source_lcao/setup_dm.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_pw/module_pwdft/VL_in_pw.h"
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/setup_exx.h"
#include "source_lcao/setup_deepks.h"

namespace LCAO_domain
{

/**
 * @brief set up wave functions, occupation numbers,
 * density matrix and charge density 
 */
template <typename TK>
void set_psi_occ_dm_chg(
		const K_Vectors &kv, // k-points
		psi::Psi<TK>* &psi, // coefficients of NAO basis
		const Parallel_Orbitals &pv, // parallel scheme of NAO basis
		elecstate::ElecState* pelec, // eigen values and weights
		LCAO_domain::Setup_DM<TK> &dmat, // density matrix 
		Charge &chr, // charge density 
		const Input_para& inp); // input parameters

/**
 * @brief set up potentials, including local pseudopotentials,
 * +U potential, solvent potential, exx potential and deepks potential 
 */
template <typename TK>
void set_pot(
        UnitCell &ucell,
		K_Vectors &kv, 
	    Structure_Factor& sf,	
		const ModulePW::PW_Basis &pw_rho, 
		const ModulePW::PW_Basis &pw_rhod, 
		elecstate::ElecState* pelec,
		const LCAO_Orbitals& orb,
		Parallel_Orbitals &pv, 
		pseudopot_cell_vl &locpp, 
        Plus_U &dftu,
        surchem& solvent,
        Exx_NAO<TK> &exx_nao,
        Setup_DeePKS<TK> &deepks,
        const Input_para &inp);

} // end namespace

#endif
