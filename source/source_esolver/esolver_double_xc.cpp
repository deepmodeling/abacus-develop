#include "esolver_double_xc.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_hamilt/module_ewald/H_Ewald_pw.h"
#include "source_hamilt/module_vdw/vdw.h"
#ifdef __MLALGO
#include "source_lcao/module_deepks/LCAO_deepks.h"
#include "source_lcao/module_deepks/LCAO_deepks_interface.h"
#endif
//-----force& stress-------------------
#include "source_lcao/FORCE_STRESS.h"

//-----HSolver ElecState Hamilt--------
#include "source_estate/elecstate_lcao.h"
#include "source_estate/elecstate_tools.h"
#include "source_lcao/hamilt_lcao.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/write_HS.h" // use ModuleIO::write_hsk()

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_DoubleXC<TK, TR>::ESolver_DoubleXC()
{
    this->classname = "ESolver_DoubleXC";
    this->basisname = "LCAO";
}

template <typename TK, typename TR>
ESolver_DoubleXC<TK, TR>::~ESolver_DoubleXC()
{
    delete this->psi_base;
    delete this->p_hamilt_base;
    delete this->pelec_base;
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "before_all_runners");
    ModuleBase::timer::tick("ESolver_DoubleXC", "before_all_runners");

    ESolver_KS_LCAO<TK, TR>::before_all_runners(ucell, inp);

    // init some items for base functional

    // 2) init ElecState
    if (this->pelec_base == nullptr)
    {
        this->pelec_base = new elecstate::ElecStateLCAO<TK>(&(this->chr_base), // use which parameter?
                                                       &(this->kv),
                                                       this->kv.get_nks(),
                                                       &(this->GG),
                                                       &(this->GK),
                                                       this->pw_rho,
                                                       this->pw_big);
    }    

    // 4) initialize electronic wave function psi
    if (this->psi_base == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = PARAM.inp.nspin;
            ncol = this->pv.ncol_bands;
            if (PARAM.inp.ks_solver == "genelpa" || PARAM.inp.ks_solver == "elpa" || PARAM.inp.ks_solver == "lapack"
                || PARAM.inp.ks_solver == "pexsi" || PARAM.inp.ks_solver == "cusolver"
                || PARAM.inp.ks_solver == "cusolvermp")
            {
                ncol = this->pv.ncol;
            }
        }
        else
        {
            nsk = this->kv.get_nks();
#ifdef __MPI
            ncol = this->pv.ncol_bands;
#else
            ncol = PARAM.inp.nbands;
#endif
        }
        this->psi_base = new psi::Psi<TK>(nsk, ncol, this->pv.nrow, this->kv.ngk, true);
    }

    // 6) initialize the density matrix
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec_base)->init_DM(&this->kv, &(this->pv), PARAM.inp.nspin);

    // 10) inititlize the charge density
    this->chr_base.allocate(PARAM.inp.nspin);
    this->pelec_base->omega = ucell.omega;

    // 11) initialize the potential
    if (this->pelec_base->pot == nullptr)
    {
        this->pelec_base->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &(this->locpp.vloc),
                                                    &(this->sf),
                                                    &(this->solvent),
                                                    &(this->pelec_base->f_en.etxc),
                                                    &(this->pelec_base->f_en.vtxc));
    }

    ModuleBase::timer::tick("ESolver_DoubleXC", "before_all_runners");
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "before_scf");
    ModuleBase::timer::tick("ESolver_DoubleXC", "before_scf");

    ESolver_KS_LCAO<TK,TR>::before_scf(ucell, istep);

    this->pelec_base->omega = ucell.omega;
    //----------------------------------------------------------
    //! calculate D2 or D3 vdW
    //----------------------------------------------------------
    auto vdw_solver = vdw::make_vdw(ucell, PARAM.inp, &(GlobalV::ofs_running));
    if (vdw_solver != nullptr)
    {
        this->pelec_base->f_en.evdw = vdw_solver->get_energy();
    }

    //----------------------------------------------------------
    //! calculate ewald energy
    //----------------------------------------------------------
    if (!PARAM.inp.test_skip_ewald)
    {
        //this->pelec_base->f_en.ewald_energy = H_Ewald_pw::compute_ewald(ucell, this->pw_rhod, this->sf.strucFac);
        this->pelec_base->f_en.ewald_energy = this->pelec->f_en.ewald_energy;
    }    

    if (this->p_hamilt_base != nullptr)
    {
        delete this->p_hamilt_base;
        this->p_hamilt_base = nullptr;
    }
    if (this->p_hamilt_base == nullptr)
    {
        elecstate::DensityMatrix<TK, double>* DM = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec_base)->get_DM();

        this->p_hamilt_base = new hamilt::HamiltLCAO<TK, TR>(
            PARAM.globalv.gamma_only_local ? &(this->GG) : nullptr,
            PARAM.globalv.gamma_only_local ? nullptr : &(this->GK),
            ucell,
            this->gd,
            &this->pv,
            this->pelec_base->pot,
            this->kv,
            this->two_center_bundle_,
            this->orb_,
            DM,
            this->deepks
#ifdef __EXX
            ,
            istep,
            GlobalC::exx_info.info_ri.real_number ? &this->exx_nao.exd->two_level_step : &this->exx_nao.exc->two_level_step,
            GlobalC::exx_info.info_ri.real_number ? &this->exx_nao.exd->get_Hexxs() : nullptr,
            GlobalC::exx_info.info_ri.real_number ? nullptr : &this->exx_nao.exc->get_Hexxs()
#endif
        );
    }

    XC_Functional::set_xc_type(PARAM.inp.deepks_out_base);
    this->pelec_base->init_scf(istep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func); 

    // DMR should be same size with Hamiltonian(R)
    dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec_base)
    ->get_DM()
    ->init_DMR(*(dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt_base)->getHR()));

    if (istep > 0)
    {
        dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec_base)->get_DM()->cal_DMR();
    }

    ModuleBase::timer::tick("ESolver_DoubleXC", "before_scf");
    return;    
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "iter_finish");
    ModuleBase::timer::tick("ESolver_DoubleXC", "iter_finish");

    bool output_iter = PARAM.inp.deepks_out_labels >0 && PARAM.inp.deepks_out_freq_elec && 
                  (iter % PARAM.inp.deepks_out_freq_elec == 0);

    if ( output_iter )
    {
        // save output charge density (density after diagnonalization)
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            ModuleBase::GlobalFunc::DCOPY(this->chr.rho[is], this->chr_base.rho[is], this->chr.rhopw->nrxx);
            if (XC_Functional::get_ked_flag())
            {
                ModuleBase::GlobalFunc::DCOPY(this->chr.kin_r[is], this->chr_base.kin_r[is], this->chr.rhopw->nrxx);
            }
        }        
    }

    ESolver_KS_LCAO<TK, TR>::iter_finish(ucell, istep, iter, conv_esolver);

    // for deepks, output labels during electronic steps (after conv_esolver is renewed)
    if ( output_iter)
    {
        // ---------- update etot and htot ----------
        // get etot of output charge density, now the etot is of density after charge mixing
        this->pelec->pot->update_from_charge(&this->chr_base, &ucell); 
        this->pelec->f_en.descf = 0.0;
        this->pelec->cal_energies(2);
        // std::cout<<"in deepks etot------"<<std::endl;
        // this->pelec->f_en.print_all();
        // std::cout<<"in deepks etot------"<<std::endl;
        // GlobalV::ofs_running << std::setprecision(15) << " in deepks etot: etot of target functional (Ry) " << this->pelec->f_en.etot << std::endl;

        // update p_hamilt using output charge density
        // Note!!!
        // This will change the result of out_mat_hs
        // The original result of out_mat_hs is H of input density, but this change H to that of output density
        // When converged, these two should be close
        if (PARAM.inp.deepks_v_delta > 0 && PARAM.inp.vl_in_h)
        {
            // update real space Hamiltonian
            this->p_hamilt->refresh();
        }

#ifdef __MLALGO
        // ---------- output tot and precalc ----------
        hamilt::HamiltLCAO<TK, TR>* p_ham_deepks = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);
        std::shared_ptr<LCAO_Deepks<TK>> ld_shared_ptr(&this->deepks.ld, [](LCAO_Deepks<TK>*) {});
        LCAO_Deepks_Interface<TK, TR> deepks_interface(ld_shared_ptr);

        deepks_interface.out_deepks_labels(this->pelec->f_en.etot,
                                            this->kv.get_nks(),
                                            ucell.nat,
                                            PARAM.globalv.nlocal,
                                            this->pelec->ekb,
                                            this->kv.kvec_d,
                                            ucell,
                                            this->orb_,
                                            this->gd,
                                            &(this->pv),
                                            *(this->psi),
                                            dynamic_cast<const elecstate::ElecStateLCAO<TK>*>(this->pelec)->get_DM(),
                                            p_ham_deepks,
                                            iter,
                                            conv_esolver,
                                            GlobalV::MY_RANK,
                                            GlobalV::ofs_running);
#endif
                                            
        // restore to density after charge mixing
        this->pelec->pot->update_from_charge(&this->chr, &ucell); 

        // ---------- prepare for base ----------
        // set as base functional Temporarily
        XC_Functional::set_xc_type(PARAM.inp.deepks_out_base);

        // update pot of pelec_base according to chr_base
        if (!conv_esolver)
        {
            this->pelec_base->pot->update_from_charge(&this->chr_base, &ucell);
        }
        else
        {
            this->pelec_base->cal_converged();
        }

        // ---------- e_base ----------
        // ebase use the same output density with etot, just different in xc
        this->pelec_base->f_en.eband = this->pelec->f_en.eband;
        this->pelec_base->f_en.deband = this->pelec->f_en.deband;
        this->pelec_base->f_en.demet = this->pelec->f_en.demet;
        this->pelec_base->f_en.descf = 0.0; // set descf to 0
        this->pelec_base->cal_energies(2); // 2 means Kohn-Sham functional
        // std::cout<<"in double_xc------"<<std::endl;
        // this->pelec_base->f_en.print_all();
        // std::cout<<"in double_xc------"<<std::endl;
        // GlobalV::ofs_running << std::setprecision(15) << " etot of base functional (Ry) " << pelec_base->f_en.etot << std::endl;

#ifdef __MLALGO        
        const std::string file_ebase = deepks_interface.get_filename("ebase", PARAM.inp.deepks_out_labels, iter);
        LCAO_deepks_io::save_npy_e(pelec_base->f_en.etot, file_ebase, GlobalV::MY_RANK);
#endif

        // ---------- h_base ----------
        if (PARAM.inp.deepks_v_delta > 0)
        {
            if (PARAM.inp.vl_in_h)
            {
                // update real space Hamiltonian
                this->p_hamilt_base->refresh();
            }

            // Note!!!
            // should not use ModuleIO::write_hsk() to output h_base, because it will call get_hs_pointers()
            // which will change the hsolver::DiagoElpa<double>::DecomposedState, influencing the following SCF steps     

#ifdef __MLALGO
            using TH = std::conditional_t<std::is_same<TK, double>::value, ModuleBase::matrix, ModuleBase::ComplexMatrix>;
            hamilt::HamiltLCAO<TK, TR>* p_ham_deepks_base = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt_base);
            int nks = this->kv.get_nks();
            std::vector<TH> h_tot(nks);
            DeePKS_domain::get_h_tot<TK, TH, TR>(this->pv, p_ham_deepks_base, h_tot, PARAM.globalv.nlocal, nks, 'H');

            const std::string file_htot = deepks_interface.get_filename("hbase", PARAM.inp.deepks_out_labels, iter);
            LCAO_deepks_io::save_npy_h<TK, TH>(h_tot, file_htot, PARAM.globalv.nlocal, nks, GlobalV::MY_RANK);
#endif
        }

        // ---------- o_base ----------
        if ( PARAM.inp.deepks_bandgap > 0 )
        {
            // obase isn't implemented yet
            // don't need to solve p_hamilt_base
            // just dm*p_hamilt_base, similar to cal_o_delta           
        }
    
        // restore to original xc
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func); 

    }
    // ---------- prepare for f_base ----------
    else if ( PARAM.inp.cal_force && conv_esolver )
    {
        // vnew must be updated for force_scc() even if not output_iter
        // set as base functional Temporarily
        XC_Functional::set_xc_type(PARAM.inp.deepks_out_base);
        this->pelec_base->cal_converged();
        // restore to original xc
        XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func); 
    }
    
    if ( PARAM.inp.cal_force )
    {
        if ( ! conv_esolver )
        {
            // use chr after mixing to restore veff, useful for vnew when converged
            this->pelec_base->pot->update_from_charge(&this->chr, &ucell); 
        }
        else
        {
            // copy charge
            for (int is = 0; is < PARAM.inp.nspin; is++)
            {
                ModuleBase::GlobalFunc::DCOPY(this->chr.rho[is], this->chr_base.rho[is], this->chr.rhopw->nrxx);
                if (XC_Functional::get_ked_flag())
                {
                    ModuleBase::GlobalFunc::DCOPY(this->chr.kin_r[is], this->chr_base.kin_r[is], this->chr.rhopw->nrxx);
                }
            }

            // copy dm
            int nks = this->kv.get_nks();
            auto _pes_lcao_base = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec_base);
            auto _pes_lcao = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
            for (int ik = 0; ik < nks; ik++)
            {
                _pes_lcao_base->get_DM()->set_DMK_pointer(ik, _pes_lcao->get_DM()->get_DMK_pointer(ik));
            }
            _pes_lcao_base->get_DM()->cal_DMR();
            _pes_lcao_base->ekb = _pes_lcao->ekb;
            _pes_lcao_base->wg = _pes_lcao->wg;          
        }        
    }
    ModuleBase::timer::tick("ESolver_DoubleXC", "iter_finish");
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::cal_force(UnitCell& ucell, ModuleBase::matrix& force)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "cal_force");
    ModuleBase::timer::tick("ESolver_DoubleXC", "cal_force");

    ModuleBase::matrix force_base;
    ModuleBase::matrix stress_base;

    Force_Stress_LCAO<TK> fsl(this->RA, ucell.nat);

    // set as base functional Temporarily
    XC_Functional::set_xc_type(PARAM.inp.deepks_out_base);

    fsl.getForceStress(ucell,
                       PARAM.inp.cal_force,
                       PARAM.inp.cal_stress,
                       PARAM.inp.test_force,
                       PARAM.inp.test_stress,
                       this->gd,
                       this->pv,
                       this->pelec_base,
                       this->psi,
                       this->GG, // mohan add 2024-04-01
                       this->GK, // mohan add 2024-04-01
                       this->two_center_bundle_,
                       this->orb_,
                       force_base,
                       stress_base,
                       this->locpp,
                       this->sf,
                       this->kv,
                       this->pw_rho,
                       this->solvent,
#ifdef __MLALGO
                       this->deepks.ld,
                       "base",
#endif
					   this->exx_nao,
					   &ucell.symm);
    // restore to original xc
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func); 

    // this will delete RA, so call it later
    ESolver_KS_LCAO<TK, TR>::cal_force(ucell, force);

    ModuleBase::timer::tick("ESolver_DoubleXC", "cal_force");
}

template class ESolver_DoubleXC<double, double>;
template class ESolver_DoubleXC<std::complex<double>, double>;
template class ESolver_DoubleXC<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver
