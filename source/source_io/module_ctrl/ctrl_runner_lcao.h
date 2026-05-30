#ifndef CTRL_RUNNER_LCAO_H 
#define CTRL_RUNNER_LCAO_H 

#include "source_cell/unitcell.h" // use UnitCell
#include "source_cell/klist.h" // use K_Vectors
#include "source_estate/elecstate.h" // use elecstate::ElecStateLCAO<TK> 
#include "source_psi/psi.h" // use Psi<TK>
#include "source_lcao/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>
#include "source_basis/module_nao/two_center_bundle.h" // use TwoCenterBundle
#include "source_lcao/setup_exx.h" // for exx, mohan add 20251018
#include "source_lcao/setup_dm.h" // for density matrix, mohan add 20251103

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_runner_lcao(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecState* pelec,// electronic info
        const LCAO_domain::Setup_DM<TK> &dmat, // mohan add 2025-11-02
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<TK>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<TK, TR>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
		ModuleBase::matrix &vloc,     // local pseudopotential 
		Exx_NAO<TK> &exx_nao,
		surchem &solvent);             // solvent model

}

#endif
