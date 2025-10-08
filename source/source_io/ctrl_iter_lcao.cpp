#include "source_io/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 
#include "source_pw/module_pwdft/global.h" // use GlobalC::restart

#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_deepks/LCAO_deepks_interface.h"
#endif

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_iter_lcao(UnitCell& ucell, // unit cell *
        const Input_para& inp, // input parameters *
		K_Vectors& kv, // k points *
		elecstate::ElecStateLCAO<TK>* pelec, // electronic info * 
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<TK>* psi, // wave functions *
        Charge &chr, // charge density *
        Charge_Mixing* p_chgmix, // charge mixing *
		hamilt::HamiltLCAO<TK, TR>* p_hamilt, // hamiltonian *
		LCAO_Orbitals &orb, // orbital info *
#ifdef __MLALGO
		LCAO_Deepks<TK>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<TK, double>& exd,
		Exx_LRI_Interface<TK, std::complex<double>>& exc,
#endif
        int &iter,
        const int istep,
        bool &conv_esolver,
		const double &scf_ene_thr)
{
    ModuleBase::TITLE("ModuleIO", "ctrl_iter_lcao");
    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_lcao");

    // save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < inp.nspin; ++is)
        {
            GlobalC::restart.save_disk("charge", is, chr.nrxx, chr.rho[is]);
        }
    }

#ifdef __EXX
    // save exx matrix
    if (inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            GlobalC::exx_info.info_ri.real_number ?
              exd.exx_iter_finish(kv, ucell, *p_hamilt, *pelec,
                *p_chgmix, scf_ene_thr, iter, istep, conv_esolver) :
              exc.exx_iter_finish(kv, ucell, *p_hamilt, *pelec,
                *p_chgmix, scf_ene_thr, iter, istep, conv_esolver);
        }
    }
#endif

    // for deepks, output labels during electronic steps (after conv_esolver is renewed)
#ifdef __MLALGO
    if (inp.deepks_out_labels >0 && inp.deepks_out_freq_elec)
    {
        if (iter % inp.deepks_out_freq_elec == 0 )
        {
            std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&ld, [](LCAO_Deepks<TK>*) {});
            LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

            deepks_interface.out_deepks_labels(pelec->f_en.etot, kv.get_nks(),
              ucell.nat, PARAM.globalv.nlocal, pelec->ekb, kv.kvec_d,
              ucell, orb, gd, &pv, *psi, pelec->get_DM(),
              p_hamilt, iter, conv_esolver, GlobalV::MY_RANK, GlobalV::ofs_running);
        }
    }
#endif

    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_lcao");
}

// TK: double  TR: double 
template void ctrl_iter_lcao<double, double>(UnitCell& ucell, // unit cell *
        const Input_para& inp, // input parameters *
		K_Vectors& kv, // k points *
		elecstate::ElecStateLCAO<double>* pelec, // electronic info * 
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<double>* psi, // wave functions *
        Charge &chr, // charge density *
        Charge_Mixing* p_chgmix, // charge mixing *
		hamilt::HamiltLCAO<double, double>* p_hamilt, // hamiltonian *
		LCAO_Orbitals &orb, // orbital info *
#ifdef __MLALGO
		LCAO_Deepks<double>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<double, double>& exd,
		Exx_LRI_Interface<double, std::complex<double>>& exc,
#endif
        int &iter,
        const int istep,
        bool &conv_esolver,
		const double &scf_ene_thr);

// TK: complex<double>  TR: double 
template void ctrl_iter_lcao<std::complex<double>, double>(UnitCell& ucell, // unit cell *
        const Input_para& inp, // input parameters *
		K_Vectors& kv, // k points *
		elecstate::ElecStateLCAO<std::complex<double>>* pelec, // electronic info * 
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<std::complex<double>>* psi, // wave functions *
        Charge &chr, // charge density *
        Charge_Mixing* p_chgmix, // charge mixing *
		hamilt::HamiltLCAO<std::complex<double>, double>* p_hamilt, // hamiltonian *
		LCAO_Orbitals &orb, // orbital info *
#ifdef __MLALGO
		LCAO_Deepks<std::complex<double>>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<std::complex<double>, double>& exd,
		Exx_LRI_Interface<std::complex<double>, std::complex<double>>& exc,
#endif
        int &iter,
        const int istep,
        bool &conv_esolver,
		const double &scf_ene_thr);

// TK: complex<double>  TR: complex<double> 
template void ctrl_iter_lcao<std::complex<double>, std::complex<double>>(UnitCell& ucell, // unit cell *
        const Input_para& inp, // input parameters *
		K_Vectors& kv, // k points *
		elecstate::ElecStateLCAO<std::complex<double>>* pelec, // electronic info * 
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<std::complex<double>>* psi, // wave functions *
        Charge &chr, // charge density *
        Charge_Mixing* p_chgmix, // charge mixing *
		hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_hamilt, // hamiltonian *
		LCAO_Orbitals &orb, // orbital info *
#ifdef __MLALGO
		LCAO_Deepks<std::complex<double>>& ld,
#endif
#ifdef __EXX
		Exx_LRI_Interface<std::complex<double>, double>& exd,
		Exx_LRI_Interface<std::complex<double>, std::complex<double>>& exc,
#endif
        int &iter,
        const int istep,
        bool &conv_esolver,
		const double &scf_ene_thr);

} // end ModuleIO
