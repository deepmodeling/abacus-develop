#ifndef CTRL_ITER_LCAO_H 
#define CTRL_ITER_LCAO_H 

#include "source_cell/unitcell.h" // use UnitCell
#include "source_cell/klist.h" // use K_Vectors
#include "source_estate/elecstate_lcao.h" // use elecstate::ElecStateLCAO<TK> 
#include "source_psi/psi.h" // use Psi<TK>
#include "source_estate/module_charge/charge.h" // use charge
#include "source_estate/module_charge/charge_mixing.h" // use charge mixing
#include "source_lcao/hamilt_lcao.h" // use hamilt::HamiltLCAO<TK, TR>
#include "source_lcao/setup_exx.h" // mohan add 20251008
#include "source_lcao/setup_deepks.h" // mohan add 20251010

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_iter_lcao(UnitCell& ucell, // unit cell *
        const Input_para& inp, // input parameters *
		K_Vectors& kv, // k points *
		elecstate::ElecState* pelec, // electronic info * 
        elecstate::DensityMatrix<TK, double>& dm, // density matrix, mohan add 2025-11-03
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<TK>* psi, // wave functions *
        Charge &chr, // charge density *
        Charge_Mixing* p_chgmix, // charge mixing *
		hamilt::HamiltLCAO<TK, TR>* p_hamilt, // hamiltonian *
		LCAO_Orbitals &orb, // orbital info *
        Setup_DeePKS<TK> &deepks,
        Exx_NAO<TK> &exx_nao,
        int &iter,
        const int istep,
        bool &conv_esolver,
		const double &scf_ene_thr);

}
#endif
