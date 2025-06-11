#ifndef CTRL_OUTPUT_LCAO_H 
#define CTRL_OUTPUT_LCAO_H 

#include <complex>

#include "module_cell/unitcell.h" // use UnitCell
#include "module_cell/klist.h" // use K_Vectors
#include "module_elecstate/elecstate_lcao.h" // use elecstate::ElecStateLCAO<TK> 
#include "module_psi/psi.h" // use Psi<TK>
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>
#include "module_basis/module_nao/two_center_bundle.h" // use TwoCenterBundle
#include "module_hamilt_lcao/module_gint/gint_k.h" // use Gint_k
#include "module_basis/module_pw/pw_basis_k.h" // use ModulePW::PW_Basis_K and ModulePW::PW_Basis
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h" // use Structure_Factor 
#include "module_rdmft/rdmft.h" // use RDMFT codes
#ifdef __EXX
#include "module_ri/Exx_LRI_interface.h" // use EXX codes
#endif

namespace ModuleIO
{
    // in principle, we need to add const for all of the variables, mohan note 2025-06-05
	template <typename TK, typename TR>
		void ctrl_output_lcao(UnitCell& ucell, 
				K_Vectors& kv,
				elecstate::ElecStateLCAO<TK>* pelec, 
				Parallel_Orbitals& pv,
				Grid_Driver& gd,
				psi::Psi<TK>* psi,
				hamilt::HamiltLCAO<TK, TR>* p_hamilt,
				TwoCenterBundle &two_center_bundle,
				Gint_k &gk,
				LCAO_Orbitals &orb,
				const ModulePW::PW_Basis_K* pw_wfc, // for berryphase
				const ModulePW::PW_Basis* pw_rho, // for berryphase
				Grid_Technique &gt, // for berryphase
				const ModulePW::PW_Basis_Big* pw_big, // for Wannier90
				const Structure_Factor& sf, // for Wannier90
				rdmft::RDMFT<TK, TR> &rdmft_solver, // for RDMFT
#ifdef __DEEPKS
				LCAO_Deepks<TK>& ld,
#endif
#ifdef __EXX
				Exx_LRI_Interface<TK, double>& exd,
				Exx_LRI_Interface<TK, std::complex<double>>& exc,
#endif
				const int istep);
}
#endif
