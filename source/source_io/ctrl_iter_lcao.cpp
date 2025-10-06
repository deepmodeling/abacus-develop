#include "source_io/ctrl_iter_lcao.h" // use ctrl_iter_lcao() 

namespace ModuleIO
{

template <typename TK, typename TR>
void ctrl_iter_lcao(UnitCell& ucell,
        const Input_para& inp,
		K_Vectors& kv, // k points *
		elecstate::ElecStateLCAO<TK>* pelec, // electronic info * 
		Parallel_Orbitals& pv, // parallel orbital info *
		Grid_Driver& gd, // adjacent atom info *
		psi::Psi<TK>* psi, // wave functions *
		hamilt::HamiltLCAO<TK, TR>* p_hamilt,
		LCAO_Orbitals &orb, // orbital info *
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

    // save charge density
    // Peize Lin add 2020.04.04
    if (GlobalC::restart.info_save.save_charge)
    {
        for (int is = 0; is < PARAM.inp.nspin; ++is)
        {
            GlobalC::restart.save_disk("charge", is, this->chr.nrxx, this->chr.rho[is]);
        }
    }

#ifdef __EXX
    // save exx matrix
    if (PARAM.inp.calculation != "nscf")
    {
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            GlobalC::exx_info.info_ri.real_number ?
              exd->exx_iter_finish(kv, ucell, *this->p_hamilt, *pelec,
                *this->p_chgmix, this->scf_ene_thr, iter, istep, conv_esolver) :
              exc->exx_iter_finish(kv, ucell, *this->p_hamilt, *pelec,
                *this->p_chgmix, this->scf_ene_thr, iter, istep, conv_esolver);
        }
    }
#endif

    // for deepks, output labels during electronic steps (after conv_esolver is renewed)
#ifdef __MLALGO
    if (PARAM.inp.deepks_out_labels >0 && PARAM.inp.deepks_out_freq_elec)
    {
        if (iter % PARAM.inp.deepks_out_freq_elec == 0 )
        {
            hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
            std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&ld, [](LCAO_Deepks<TK>*) {});
            LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

            deepks_interface.out_deepks_labels(pelec->f_en.etot, kv.get_nks(),
              ucell.nat, PARAM.globalv.nlocal, pelec->ekb, kv.kvec_d,
              ucell, orb, gd, &pv, *psi,
              dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(pelec)->get_DM(),
              p_ham_deepks, iter, conv_esolver, GlobalV::MY_RANK, GlobalV::ofs_running);
        }
    }
#endif

    ModuleBase::timer::tick("ModuleIO", "ctrl_iter_lcao");
}
