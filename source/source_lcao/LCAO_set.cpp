#include "source_lcao/LCAO_set.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/setup_psi.h" // use Setup_Psi
#include "source_io/module_wf/read_wfc_nao.h" // use read_wfc_nao
#include "source_estate/elecstate_tools.h" // use fixed_weights
#include "source_lcao/module_hcontainer/read_hcontainer.h"
#include "source_lcao/rho_tau_lcao.h" // use dm2rho
#include "source_lcao/hamilt_lcao.h" // use HamiltLCAO for init_chg_hr
#include "source_hsolver/hsolver_lcao.h" // use HSolverLCAO for init_chg_hr

template <typename TK>
void LCAO_domain::set_psi_occ_dm_chg(
		const K_Vectors &kv, // k-points
		psi::Psi<TK>* &psi, // coefficients of NAO basis
		const Parallel_Orbitals &pv, // parallel scheme of NAO basis
		elecstate::ElecState* pelec, // eigen values and weights
		LCAO_domain::Setup_DM<TK> &dmat, // density matrix 
		Charge &chr, // charge density 
		const Input_para &inp) // input parameters
{

    //! 1) init electronic wave function psi
    Setup_Psi<TK>::allocate_psi(psi, kv, pv, inp);

    //! 2) read psi from file
    if (inp.init_wfc == "file" && inp.esolver_type != "tddft")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
             pv, *psi, pelec->ekb, pelec->wg, kv.ik2iktot,
             kv.get_nkstot(), inp.nspin))
        {
            ModuleBase::WARNING_QUIT("set_psi_occ_dm_chg", "read electronic wave functions failed");
        }
    }

    //! 3) set occupations, tddft does not need to set occupations in the first scf
    if (inp.ocp && inp.esolver_type != "tddft")
    {
        elecstate::fixed_weights(inp.ocp_kb, inp.nbands, inp.nelec,
          &kv, pelec->wg, pelec->skip_weights);
    }

    //! 4) init DMK, but DMR is constructed in before_scf()
    dmat.allocate_dm(&kv, &pv, inp.nspin);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "CHARGE");

    return;
}


template <typename TK>
void LCAO_domain::set_pot(
        UnitCell &ucell, // not const because of dftu
		K_Vectors &kv, // not const due to exx 
	    Structure_Factor& sf, // will be modified in potential	
		const ModulePW::PW_Basis &pw_rho, 
		const ModulePW::PW_Basis &pw_rhod, 
		elecstate::ElecState* pelec,
		const LCAO_Orbitals& orb,
		Parallel_Orbitals &pv, // not const due to deepks 
		pseudopot_cell_vl &locpp, 
        Plus_U &dftu,
        surchem& solvent,
        Exx_NAO<TK> &exx_nao,
        Setup_DeePKS<TK> &deepks,
        const Input_para &inp)
{
    //! 1) init local pseudopotentials
    locpp.init_vloc(ucell, &pw_rho);

    //! 2) init potentials
    if (pelec->pot == nullptr)
    {
        // where is the pot deleted?
        pelec->pot = new elecstate::Potential(&pw_rhod, &pw_rho,
          &ucell, &locpp.vloc, &sf, &solvent,
          &(pelec->f_en.etxc), &(pelec->f_en.vtxc));
    }

    //! 3) initialize DFT+U
    if (inp.dft_plus_u)
    {
        dftu.init(ucell, &pv, kv.get_nks(), &orb);
    }

    //! 4) init exact exchange calculations
    exx_nao.before_runner(ucell, kv, orb, pv, inp);

    //! 5) init deepks
    deepks.before_runner(ucell, kv.get_nks(), orb, pv, inp);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "POTENTIALS");

    return;
}

template <typename TK>
void LCAO_domain::init_dm_from_file(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<TK>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv)
{
    ModuleBase::TITLE("LCAO_domain", "init_dm_from_file");
    const int nspin_dm = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin_dm; ++is)
    {
        const std::string dmfile = readin_dir + "/dmrs" + std::to_string(is + 1) + "_nao.csr";
        hamilt::HContainer<double>* dm_container = dmat.dm->get_DMR_vector()[is];
        hamilt::Read_HContainer<double> reader_dm(
            dm_container,
            dmfile,
            PARAM.globalv.nlocal,
            &ucell
        );
        reader_dm.read();
    }
    return;
}

template <typename TK>
void LCAO_domain::init_chg_dm(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<TK>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    Charge* chr)
{
    ModuleBase::TITLE("LCAO_domain", "init_chg_dm");

    // Step 1: Read density matrix from file
    LCAO_domain::init_dm_from_file<TK>(readin_dir, nspin, dmat, ucell, pv);

    // Step 2: Convert density matrix to charge density
    LCAO_domain::dm2rho(dmat.dm->get_DMR_vector(), nspin, chr, true);

    return;
}

template <typename TR>
void LCAO_domain::init_hr_from_file(
    const std::string hrfile,
    hamilt::HContainer<TR>* hmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv)
{
    ModuleBase::TITLE("LCAO_domain", "init_hr_from_file");

    // Check if file exists
    std::ifstream test_file(hrfile);
    if (!test_file.good())
    {
        std::string error_msg = "Cannot open Hamiltonian file: " + hrfile + "\n\n";
        error_msg += "When using init_chg=hr, you need to provide Hamiltonian matrix files:\n";
        error_msg += "  - For nspin=1: hrs1_nao.csr\n";
        error_msg += "  - For nspin=2: hrs1_nao.csr (spin-up) and hrs2_nao.csr (spin-down)\n\n";
        error_msg += "Solutions:\n";
        error_msg += "  1. Run an SCF calculation first with 'out_mat_hs2 1' to generate HR files\n";
        error_msg += "  2. Check that 'read_file_dir' points to the correct directory\n";
        error_msg += "  3. Use 'init_chg file' or 'init_chg atomic' instead";
        ModuleBase::WARNING_QUIT("LCAO_domain::init_hr_from_file", error_msg);
    }
    test_file.close();

    hmat->set_zero();
    hamilt::Read_HContainer<TR> reader_hr(hmat, hrfile, PARAM.globalv.nlocal, &ucell);
    reader_hr.read();
    return;
}

template <typename TK, typename TR>
void LCAO_domain::init_chg_hr(
    const std::string& readin_dir,
    const int nspin,
    hamilt::Hamilt<TK>* p_hamilt,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    psi::Psi<TK>& psi,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<TK, double>& dm,
    Charge& chr,
    const std::string& ks_solver)
{
    ModuleBase::TITLE("LCAO_domain", "init_chg_hr");

    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(p_hamilt);
    if (!hamilt_lcao)
    {
        ModuleBase::WARNING_QUIT("LCAO_domain::init_chg_hr", "p_hamilt is not HamiltLCAO");
    }

    // Step 1: Read HR from file(s)
    if (nspin == 2)
    {
        // nspin=2: load spin-up into first half of hRS2, spin-down into second half
        const std::string hrfile_up = readin_dir + "/hrs1_nao.csr";
        LCAO_domain::init_hr_from_file<TR>(hrfile_up, hamilt_lcao->getHR(), ucell, pv);

        // switch hR data pointer to spin-down half, then read hrs2
        auto& hRS2 = hamilt_lcao->getHRS2();
        hamilt_lcao->getHR()->allocate(hRS2.data() + hRS2.size() / 2, 0);
        const std::string hrfile_down = readin_dir + "/hrs2_nao.csr";
        LCAO_domain::init_hr_from_file<TR>(hrfile_down, hamilt_lcao->getHR(), ucell, pv);

        // restore hR to spin-up half
        hamilt_lcao->getHR()->allocate(hRS2.data(), 0);
    }
    else
    {
        const std::string hrfile = readin_dir + "/hrs1_nao.csr";
        LCAO_domain::init_hr_from_file<TR>(hrfile, hamilt_lcao->getHR(), ucell, pv);
    }

    // Step 2: Mark HR as loaded from file (skip operator recalculation)
    p_hamilt->refresh(false);

    // Step 3: Diagonalize to get wavefunctions and charge density
    hsolver::HSolverLCAO<TK> hsolver_lcao_obj(pv, ks_solver);
    hsolver_lcao_obj.solve(p_hamilt, psi, pelec, dm, chr, nspin, 0);
}



template void LCAO_domain::set_psi_occ_dm_chg<double>(
		const K_Vectors &kv, // k-points
		psi::Psi<double>* &psi, // coefficients of NAO basis
		const Parallel_Orbitals &pv, // parallel scheme of NAO basis
		elecstate::ElecState* pelec, // eigen values and weights
		LCAO_domain::Setup_DM<double> &dmat, // density matrix 
		Charge &chr, // charge density 
		const Input_para &inp);

template void LCAO_domain::set_psi_occ_dm_chg<std::complex<double>>(
		const K_Vectors &kv, // k-points
		psi::Psi<std::complex<double>>* &psi, // coefficients of NAO basis
		const Parallel_Orbitals &pv, // parallel scheme of NAO basis
		elecstate::ElecState* pelec, // eigen values and weights
		LCAO_domain::Setup_DM<std::complex<double>> &dmat, // density matrix 
		Charge &chr, // charge density 
		const Input_para &inp);

template void LCAO_domain::set_pot<double>(
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
        Exx_NAO<double> &exx_nao,
        Setup_DeePKS<double> &deepks,
        const Input_para &inp);

template void LCAO_domain::set_pot<std::complex<double>>(
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
        Exx_NAO<std::complex<double>> &exx_nao,
        Setup_DeePKS<std::complex<double>> &deepks,
        const Input_para &inp);

template void LCAO_domain::init_dm_from_file<double>(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<double>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv);
template void LCAO_domain::init_dm_from_file<std::complex<double>>(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<std::complex<double>>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv);

template void LCAO_domain::init_chg_dm<double>(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<double>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    Charge* chr);
template void LCAO_domain::init_chg_dm<std::complex<double>>(
    const std::string& readin_dir,
    const int nspin,
    LCAO_domain::Setup_DM<std::complex<double>>& dmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    Charge* chr);

template void LCAO_domain::init_hr_from_file<double>(
    const std::string hrfile,
    hamilt::HContainer<double>* hmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv);
template void LCAO_domain::init_hr_from_file<std::complex<double>>(
    const std::string hrfile,
    hamilt::HContainer<std::complex<double>>* hmat,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv);

template void LCAO_domain::init_chg_hr<double, double>(
    const std::string& readin_dir,
    const int nspin,
    hamilt::Hamilt<double>* p_hamilt,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    psi::Psi<double>& psi,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<double, double>& dm,
    Charge& chr,
    const std::string& ks_solver);
template void LCAO_domain::init_chg_hr<std::complex<double>, double>(
    const std::string& readin_dir,
    const int nspin,
    hamilt::Hamilt<std::complex<double>>* p_hamilt,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    psi::Psi<std::complex<double>>& psi,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<std::complex<double>, double>& dm,
    Charge& chr,
    const std::string& ks_solver);
template void LCAO_domain::init_chg_hr<std::complex<double>, std::complex<double>>(
    const std::string& readin_dir,
    const int nspin,
    hamilt::Hamilt<std::complex<double>>* p_hamilt,
    const UnitCell& ucell,
    const Parallel_Orbitals* pv,
    psi::Psi<std::complex<double>>& psi,
    elecstate::ElecState* pelec,
    elecstate::DensityMatrix<std::complex<double>, double>& dm,
    Charge& chr,
    const std::string& ks_solver);
