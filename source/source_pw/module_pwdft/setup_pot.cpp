#include "source_pw/module_pwdft/setup_pot.h"

#include "source_estate/module_charge/symmetry_rho.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/onsite_projector.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_pw/module_pwdft/VSep_in_pw.h"

template <typename T, typename Device>
void pw::setup_pot(const int istep,
		UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<T, Device>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<T, Device>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp) // input parameters
{
    ModuleBase::TITLE("pw", "setup_pot");

    //----------------------------------------------------------
    //! 0) DFT-1/2 calculations, sep potential need to generate
    // before effective potential calculation
    //----------------------------------------------------------
    if (PARAM.inp.dfthalf_type > 0)
    {
        vsep_cell->generate_vsep_r(pw_rhod[0], sf.strucFac, ucell.sep_cell);
    }

    //----------------------------------------------------------
    //! 1) Renew local pseudopotential
    //----------------------------------------------------------
    pelec->init_scf(istep,
            ucell,
            para_grid,
            sf.strucFac,
            locpp.numeric,
			ucell.symm,
			(void*)pw_wfc);

	//----------------------------------------------------------
	//! 2) Symmetrize the charge density (rho)
	//----------------------------------------------------------

	//! Symmetry_rho should behind init_scf, because charge should be
	//! initialized first. liuyu comment: Symmetry_rho should be
	//! located between init_rho and v_of_rho?
	Symmetry_rho srho;
	for (int is = 0; is < inp.nspin; is++)
	{
		srho.begin(is, chr, pw_rhod, ucell.symm);
	}

	//----------------------------------------------------------
	//! 3) Calculate the effective potential with rho
	//----------------------------------------------------------
	//! liuyu move here 2023-10-09
	//! D in uspp need vloc, thus behind init_scf()
	//! calculate the effective coefficient matrix
	//! for non-local pseudopotential projectors
	ModuleBase::matrix veff = pelec->pot->get_effective_v();

	ppcell.cal_effective_D(veff, pw_rhod, ucell);

	//----------------------------------------------------------
	//! 4) Onsite projectors
	//----------------------------------------------------------
	if (PARAM.inp.onsite_radius > 0)
	{
		auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
		onsite_p->init(PARAM.inp.orbital_dir,
				&ucell,
				*(kspw_psi),
				kv,
				*(pw_wfc),
				sf,
				PARAM.inp.onsite_radius,
				PARAM.globalv.nqx,
				PARAM.globalv.dq,
				pelec->wg,
				pelec->ekb);
	}

    //----------------------------------------------------------
    //! 5) Spin-constrained algorithms
    //----------------------------------------------------------
    if (PARAM.inp.sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        sc.init_sc(PARAM.inp.sc_thr,
                   PARAM.inp.nsc,
                   PARAM.inp.nsc_min,
                   PARAM.inp.alpha_trial,
                   PARAM.inp.sccut,
                   PARAM.inp.sc_drop_thr,
                   ucell,
                   nullptr, // parallel orbitals
                   PARAM.inp.nspin,
                   kv,
                   p_hamilt,
                   kspw_psi,
#ifdef __LCAO
                   nullptr, // density matrix, not useful in LCAO, mohan note 2025-11-03
#endif
                   pelec,
                   pw_wfc);
    }

    //----------------------------------------------------------
    //! 6) DFT+U algorithm
    //----------------------------------------------------------
    if (PARAM.inp.dft_plus_u)
    {
        auto* dftu = ModuleDFTU::DFTU::get_instance();
        dftu->init(ucell, nullptr, kv.get_nks());
    }

    return;
}

template void pw::setup_pot<std::complex<float>, base_device::DEVICE_CPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<float>, base_device::DEVICE_CPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_CPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters


template void pw::setup_pot<std::complex<double>, base_device::DEVICE_CPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_CPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

#if ((defined __CUDA) || (defined __ROCM))

template void pw::setup_pot<std::complex<float>, base_device::DEVICE_GPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<float>, base_device::DEVICE_GPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<float>, base_device::DEVICE_GPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

template void pw::setup_pot<std::complex<double>, base_device::DEVICE_GPU>(
        const int istep,  // ionic step
		UnitCell& ucell, // unitcell
		const K_Vectors &kv, // kpoints
        Structure_Factor &sf, // structure factors
		elecstate::ElecState *pelec, // pointer of electrons
		const Parallel_Grid &para_grid, // parallel of FFT grids
		const Charge &chr, // charge density
		pseudopot_cell_vl &locpp, // local pseudopotentials
		pseudopot_cell_vnl &ppcell, // non-local pseudopotentials
		VSep* vsep_cell, // U-1/2 method
		psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* kspw_psi, // electronic wave functions
        hamilt::Hamilt<std::complex<double>, base_device::DEVICE_GPU>* p_hamilt, // hamiltonian
		ModulePW::PW_Basis_K *pw_wfc,  // pw for wfc
		const ModulePW::PW_Basis *pw_rhod, // pw for rhod
		const Input_para& inp); // input parameters

#endif
