#include <fstream>
#include <algorithm>

#include "esolver_sdft_pw.h"
#include "module_base/timer.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_elecstate/elecstate_pw_sdft.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_io/rho_io.h"
#include "module_io/write_istate_info.h"

//-------------------Temporary------------------
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
//----------------------------------------------
//-----force-------------------
#include "module_hamilt_pw/hamilt_stodft/sto_forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_stodft/sto_stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

ESolver_SDFT_PW::ESolver_SDFT_PW()
{
    classname = "ESolver_SDFT_PW";
    basisname = "PW";
}

ESolver_SDFT_PW::~ESolver_SDFT_PW()
{
}

void ESolver_SDFT_PW::Init(Input &inp, UnitCell &ucell)
{
    this->nche_sto = inp.nche_sto;
    ESolver_KS::Init(inp,ucell);

    
    this->pelec = new elecstate::ElecStatePW_SDFT( GlobalC::wfcpw, &(chr), (K_Vectors*)(&(GlobalC::kv)));

    // Inititlize the charge density.
    this->pelec->charge->allocate(GlobalV::NSPIN, GlobalC::rhopw->nrxx, GlobalC::rhopw->npw);

    // Initializee the potential.
    if(this->pelec->pot == nullptr)
    {
        this->pelec->pot = new elecstate::Potential(
            GlobalC::rhopw,
            &GlobalC::ucell,
            &(GlobalC::ppcell.vloc),
            &(GlobalC::sf.strucFac),
            &(GlobalC::en.etxc),
            &(GlobalC::en.vtxc)
        );
        GlobalTemp::veff = &(this->pelec->pot->get_effective_v());
    }

    //Maybe NSPIN=2 is not considered in this ESolver, but FYI
    //Fix pelec->wg by ocp_kb
    if(GlobalV::ocp)
    {
        this->pelec->fixed_weights(GlobalV::ocp_kb.data());
    }

    this->Init_GlobalC(inp,ucell);//temporary

	stowf.init(GlobalC::kv.nks);
	if(INPUT.nbands_sto != 0)	Init_Sto_Orbitals(this->stowf, inp.seed_sto, GlobalC::kv.nks);
	else						Init_Com_Orbitals(this->stowf, GlobalC::kv, GlobalC::wf.npwx);
	for (int ik =0 ; ik < GlobalC::kv.nks; ++ik)
    {
        this->stowf.shchi[ik].create(this->stowf.nchip[ik],GlobalC::wf.npwx,false);
        if(GlobalV::NBANDS > 0)
        {
            this->stowf.chiortho[ik].create(this->stowf.nchip[ik],GlobalC::wf.npwx,false);
        }
    }

    this->phsol = new hsolver::HSolverPW_SDFT(GlobalC::wfcpw, this->stowf, inp.method_sto);
   

}

void ESolver_SDFT_PW::beforescf(const int istep)
{
    ESolver_KS_PW::beforescf(istep);
	if(istep > 0 && INPUT.nbands_sto != 0 && INPUT.initsto_freq > 0 && istep%INPUT.initsto_freq == 0) Update_Sto_Orbitals(this->stowf, INPUT.seed_sto, GlobalC::kv.nks);
}

void ESolver_SDFT_PW::eachiterfinish(int iter)
{
	//this->pelec->print_eigenvalue(GlobalV::ofs_running);
    GlobalC::en.calculate_etot();
}
void ESolver_SDFT_PW::afterscf(const int istep)
{
    if(GlobalV::out_chg > 0)
    {
	    for(int is=0; is<GlobalV::NSPIN; is++)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
            double& ef_tmp = GlobalC::en.get_ef(is,GlobalV::TWO_EFERMI);
            ModuleIO::write_rho(
#ifdef __MPI
                GlobalC::bigpw->bz,
                GlobalC::bigpw->nbz,
                GlobalC::rhopw->nplane,
                GlobalC::rhopw->startz_current,
#endif
                pelec->charge->rho_save[is],
                is,
                GlobalV::NSPIN,
                0,
                ssc.str(),
                GlobalC::rhopw->nx,
                GlobalC::rhopw->ny,
                GlobalC::rhopw->nz,
                ef_tmp,
                &(GlobalC::ucell));
        }
    }
    if(this->conv_elec)
    {
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        GlobalV::ofs_running << " convergence has NOT been achieved!" << std::endl;
    }
}

void ESolver_SDFT_PW::hamilt2density(int istep, int iter, double ethr)
{
	// reset energy 
    this->pelec->eband  = 0.0;
    this->pelec->demet  = 0.0;
    this->pelec->ef     = 0.0;
    GlobalC::en.ef_up  = 0.0;
    GlobalC::en.ef_dw  = 0.0;
    // choose if psi should be diag in subspace
    // be careful that istep start from 0 and iter start from 1
    if(istep==0&&iter==1) 
    {
        hsolver::DiagoIterAssist<double>::need_subspace = false;
    }
    else 
    {
        hsolver::DiagoIterAssist<double>::need_subspace = true;
	}
    hsolver::DiagoIterAssist<double>::PW_DIAG_THR = ethr; 
    hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;
    this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec,this->stowf, istep, iter, GlobalV::KS_SOLVER);   
    if(GlobalV::MY_STOGROUP==0)
    {
        Symmetry_rho srho;
        for(int is=0; is < GlobalV::NSPIN; is++)
        {
            srho.begin(is, *(this->pelec->charge), GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
        }
        GlobalC::en.deband = GlobalC::en.delta_e(this->pelec);
    }
    else
    {
#ifdef __MPI
			if(ModuleSymmetry::Symmetry::symm_flag == 1)	MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    // transform energy for print
    GlobalC::en.eband = this->pelec->eband;
    GlobalC::en.demet = this->pelec->demet;
    GlobalC::en.ef = this->pelec->ef; 
}

void ESolver_SDFT_PW::cal_Energy(double& etot)
{
    etot = GlobalC::en.etot;
}

void ESolver_SDFT_PW::cal_Force(ModuleBase::matrix &force)
{
	Sto_Forces ff;
    ff.init(force, this->pelec->wg, this->psi, this->stowf, pelec->charge);
}
void ESolver_SDFT_PW::cal_Stress(ModuleBase::matrix &stress)
{
	Sto_Stress_PW ss;
    ss.cal_stress(stress, this->pelec->wg, this->psi, this->stowf, pelec->charge);
}
void ESolver_SDFT_PW::postprocess()
{

    GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    ModuleIO::write_istate_info(this->pelec->ekb,this->pelec->wg,&(GlobalC::kv),&(GlobalC::Pkpoints));

    if(this->maxniter == 0)
    {
        int iter = 1;
        int istep = 0;
        hsolver::DiagoIterAssist<double>::PW_DIAG_NMAX = GlobalV::PW_DIAG_NMAX;
        hsolver::DiagoIterAssist<double>::PW_DIAG_THR = std::max(std::min(1e-5, 0.1 * GlobalV::SCF_THR / std::max(1.0, GlobalV::nelec)),1e-12);
        hsolver::DiagoIterAssist<double>::need_subspace = false;
        this->phsol->solve(this->p_hamilt, this->psi[0], this->pelec,this->stowf,istep, iter, GlobalV::KS_SOLVER, true);
        GlobalC::en.ef = this->pelec->ef; //Temporary: Please use this->pelec->ef. GlobalC::en.ef is not recommended.
    }
    ((hsolver::HSolverPW_SDFT*)phsol)->stoiter.cleanchiallorder();//release lots of memories
    int nche_test = 0;
    if(INPUT.cal_cond)  nche_test = std::max(nche_test, INPUT.cond_nche);
    if(INPUT.out_dos)  nche_test = std::max(nche_test, INPUT.dos_nche);
    if(nche_test > 0)   check_che(nche_test);
    
    if(INPUT.cal_cond)
	{
        this->sKG(INPUT.cond_nche,INPUT.cond_fwhm,INPUT.cond_wcut,INPUT.cond_dw,INPUT.cond_dt, INPUT.cond_dtbatch);
    }
    if(INPUT.out_dos)
	{
        double emax, emin;
        if(INPUT.dos_setemax)	
            emax = INPUT.dos_emax_ev;
        else
            emax =  ((hsolver::HSolverPW_SDFT*)phsol)->stoiter.stohchi.Emax*ModuleBase::Ry_to_eV;
		if(INPUT.dos_setemin)
        	emin = INPUT.dos_emin_ev;
        else
            emin =  ((hsolver::HSolverPW_SDFT*)phsol)->stoiter.stohchi.Emin*ModuleBase::Ry_to_eV;
		if(!INPUT.dos_setemax && !INPUT.dos_setemin)
		{
			double delta=(emax-emin)*INPUT.dos_scale;
			emax=emax+delta/2.0;
			emin=emin-delta/2.0;
		}
        this->caldos(INPUT.dos_nche, INPUT.dos_sigma, emin, emax, INPUT.dos_edelta_ev, INPUT.npart_sto );
    }
}

}