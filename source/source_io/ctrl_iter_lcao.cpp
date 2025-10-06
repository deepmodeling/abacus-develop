#include "source_io/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_iter_lcao(UnitCell& ucell,
        const Input_para& inp,
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
#ifdef __MLALGO
		LCAO_Deepks<TK>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<TK, double>& exd,
		Exx_LRI_Interface<TK, std::complex<double>>& exc,
#endif
		const int istep)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_iter_lcao");
    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_lcao");



    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_lcao");
}
