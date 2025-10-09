#ifndef CTRL_RUNNER_LCAO_H 
#define CTRL_RUNNER_LCAO_H 

#include "source_cell/unitcell.h" // use UnitCell
#include "source_cell/klist.h" // use K_Vectors
#include "source_estate/elecstate_lcao.h" // use elecstate::ElecStateLCAO<TK> 
#include "source_psi/psi.h" // use Psi<TK>
#include "source_lcao/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>
#include "source_basis/module_nao/two_center_bundle.h" // use TwoCenterBundle
#include "source_lcao/module_gint/gint_k.h" // use Gint_k
#ifdef __EXX
#include "source_lcao/module_ri/Exx_LRI_interface.h" // use EXX codes
#endif

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_runner_lcao(UnitCell& ucell,      // unitcell
        const Input_para &inp,              // input
		K_Vectors &kv,                      // k-point
		elecstate::ElecStateLCAO<TK>* pelec,// electronic info
		Parallel_Orbitals &pv,              // orbital info
        Parallel_Grid &pgrid,               // grid info
		Grid_Driver &gd,                    // search for adjacent atoms
		psi::Psi<TK>* psi,                  // wave function
        Charge &chr,                  // charge density
		hamilt::HamiltLCAO<TK, TR>* p_hamilt, // hamiltonian
		TwoCenterBundle &two_center_bundle,   // use two-center integration
        Gint_Gamma &gg,                     // gint for Gamma-only
		Gint_k &gk,                         // gint for multi k-points
		LCAO_Orbitals &orb,                 // LCAO orbitals
		ModulePW::PW_Basis* pw_rho,   // charge density
		ModulePW::PW_Basis* pw_rhod,  // dense charge density 
		Structure_Factor &sf,         // structure factor
        ModuleBase::matrix &vloc,     // local pseudopotential 
#ifdef __EXX
		std::shared_ptr<Exx_LRI_Interface<TK, double>> exd,
		std::shared_ptr<Exx_LRI_Interface<TK, std::complex<double>>> exc,
#endif
        surchem &solvent);             // solvent model

}

#endif
