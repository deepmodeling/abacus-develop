#ifndef LCAO_SET_H
#define LCAO_SET_H

#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_psi/psi.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/setup_dm.h"
#include "source_pw/module_pwdft/structure_factor.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_surchem/surchem.h"
#include "source_pw/module_pwdft/vl_pw.h"
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

/**
 * @brief read in DMR from file, and save it into dmat
 * @param readin_dir directory containing dmrs*_nao.csr files
 * @param nspin number of spin components (1 or 2)
 */
template <typename TK>
void init_dm_from_file(
	const std::string& readin_dir,
	const int nspin,
	LCAO_domain::Setup_DM<TK>& dmat,
	const UnitCell& ucell,
	const Parallel_Orbitals* pv);

/**
 * @brief initialize charge density from density matrix file (init_chg=dm)
 * This function reads DMR from file and converts it to charge density
 * @param readin_dir directory containing dmrs*_nao.csr files
 * @param nspin number of spin components (1 or 2)
 * @param dmat density matrix object
 * @param ucell unit cell
 * @param pv parallel orbitals
 * @param chr charge density object
 */
template <typename TK>
void init_chg_dm(
	const std::string& readin_dir,
	const int nspin,
	LCAO_domain::Setup_DM<TK>& dmat,
	const UnitCell& ucell,
	const Parallel_Orbitals* pv,
	Charge* chr);

/**
 * @brief read in HR from file, and save it into hmat
 */
template <typename TK>
void init_hr_from_file(
	const std::string hrfile,
	hamilt::HContainer<TK>* hmat,
	const UnitCell& ucell,
	const Parallel_Orbitals* pv);

/**
 * @brief initialize charge density from Hamiltonian matrix file (init_chg=hr)
 * Reads HR from file(s), diagonalizes to get wavefunctions, then computes charge density.
 * For nspin=2, reads both hrs1_nao.csr (spin-up) and hrs2_nao.csr (spin-down)
 * into the two halves of HamiltLCAO::hRS2.
 * @tparam TK k-space type (double or complex<double>)
 * @tparam TR real-space type (double)
 * @param readin_dir directory containing hrs*_nao.csr files
 * @param nspin number of spin components
 * @param p_hamilt pointer to Hamilt base class (will be dynamic_cast to HamiltLCAO)
 * @param ucell unit cell
 * @param pv parallel orbitals
 * @param psi wave function object
 * @param pelec electronic state
 * @param dm density matrix
 * @param chr charge density
 * @param ks_solver solver method name
 */
template <typename TK, typename TR>
void init_chg_hr(
	const std::string& readin_dir,
	const int nspin,
	hamilt::Hamilt<TK>* p_hamilt,
	const UnitCell& ucell,
	const Parallel_Orbitals* pv,
	psi::Psi<TK>& psi,
	elecstate::ElecState* pelec,
	elecstate::DensityMatrix<TK, double>& dm,
	Charge& chr,
	const std::string& ks_solver);
} // end namespace

#endif
