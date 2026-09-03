#include "esolver_fp.h"

#include "source_cell/cal_ux.h"
#include "source_estate/module_charge/symm_rho.h"
#include "source_cell/read_pp_ucell.h"
#include "source_estate/param_update.h"
#include "source_hamilt/module_ewald/h_ewald_pw.h"
#include "source_hamilt/module_vdw/vdw.h"
#include "source_io/module_output/output_log.h"
#include "source_io/module_output/print_info.h"
#include "source_io/module_chgpot/rhog_io.h"
#include "source_io/module_parameter/parameter.h"

#include "source_pw/module_pwdft/setup_pwrho.h" // mohan 20251005
#include "source_hamilt/module_xc/xc_functional.h" // mohan 20251005
#include "source_io/module_ctrl/ctrl_output_fp.h"
#include "source_io/module_chgpot/write_init.h" // write_chg_init, write_pot_init

namespace ModuleESolver
{

ESolver_FP::ESolver_FP()
{
}

ESolver_FP::~ESolver_FP()
{
	//****************************************************
	// do not add any codes in this deconstructor funcion
	//****************************************************
    // mohan add 20251005
    pw::teardown_pwrho(this->pw_rho_flag, PARAM.globalv.double_grid, this->pw_rho, this->pw_rhod);

	delete this->pelec;
}

void ESolver_FP::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    this->inp_ = &inp;

    SurchemParameters surchem_parameters;
    surchem_parameters.eb_k = inp.eb_k;
    surchem_parameters.tau = inp.tau;
    surchem_parameters.sigma_k = inp.sigma_k;
    surchem_parameters.nc_k = inp.nc_k;
    this->solvent.set_parameters(surchem_parameters);

    ModuleBase::TITLE("ESolver_FP", "before_all_runners");

    //! 1) read pseudopotentials
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const int npol = PARAM.globalv.npol;
    const bool two_fermi = PARAM.globalv.two_fermi;
    auto atoms_info = unitcell::read_pseudo(GlobalV::ofs_running,
                                            ucell,
                                            this->inp_->pseudo_dir,
                                            global_out_dir,
                                            this->inp_->out_element_info,
                                            this->inp_->dft_functional,
                                            this->inp_->lspinorb,
                                            this->inp_->pseudo_rcut,
                                            this->inp_->soc_lambda,
                                            this->inp_->nspin,
                                            npol,
                                            this->inp_->basis_type,
                                            this->inp_->esolver_type,
                                            this->inp_->init_wfc,
                                            this->inp_->nbands,
                                            two_fermi,
                                            this->inp_->nelec_delta,
                                            this->inp_->smearing_method,
                                            this->inp_->ks_solver,
                                            this->inp_->bndpar,
                                            this->inp_->nelec,
                                            this->inp_->nupdown);
    elecstate::ParamUpdater::update_from_atoms_info(atoms_info);

    //! 2) setup pw_rho, pw_rhod, pw_big, sf, and read_pseudopotentials
    pw::setup_pwrho(ucell, PARAM.globalv.double_grid, this->pw_rho_flag, 
      this->pw_rho, this->pw_rhod, this->pw_big, this->classname, inp);

    //! 3) setup structure factors
    this->sf.set(this->pw_rhod, inp.nbspline);

    //! 4) init charge extrapolation
    this->CE.Init_CE(inp.nspin, ucell.nat, this->pw_rhod->nrxx, inp.chg_extrap);

    //! 5) symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
        ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running,
                             this->inp_->symmetry_prec, inp.nspin, this->inp_->calculation, cal_symm_repr);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");

    //! 7) setup k points in the Brillouin zone according to symmetry.
    const bool use_ibz = !inp.berry_phase && ModuleSymmetry::Symmetry::symm_flag != -1;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {this->inp_->kspacing[0], this->inp_->kspacing[1], this->inp_->kspacing[2]};
    const double koffset[3] = {this->inp_->koffset[0], this->inp_->koffset[1], this->inp_->koffset[2]};
    this->kv.set(ucell, ucell.symm, inp.kpoint_file, inp.nspin, ucell.G, ucell.latvec, GlobalV::ofs_running, use_ibz, global_out_dir, gamma_only_local, kspacing, this->inp_->kmesh_type, koffset);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");

    //! 8) print information
    ModuleIO::print_parameters(ucell, this->kv, inp);

    //! 9) parallel of FFT grid
    const int nprocgroup = (this->inp_->esolver_type == "sdft") ? GlobalV::NPROC_IN_BNDGROUP : GlobalV::NPROC;
    this->Pgrid.init(this->pw_rhod->nx, this->pw_rhod->ny, this->pw_rhod->nz,
            this->pw_rhod->nplane, this->pw_rhod->nrxx, pw_big->nbz, pw_big->bz,
            nprocgroup);

    //! 10) calculate the structure factor
    this->sf.setup(&ucell, Pgrid, this->pw_rhod);

    //! 11) setup the xc functional
    XC_Functional::set_xc_type(ucell.atoms[0].ncpp.xc_func);
    GlobalV::ofs_running<<XC_Functional::output_info()<<std::endl;

    //! 11) initialize the charge density, we need to first set xc_type,
    // then we can call chr.allocate()
	this->chr.set_rhopw(this->pw_rhod); // mohan add 20251130
    const bool kin_den = this->chr.kin_density(); // mohan add 20251202
	this->chr.allocate(inp.nspin, kin_den); // mohan move this from setup_estate_pw, 20251128


    return;
}

void ESolver_FP::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_FP", "after_scf");

    //! Output convergence information
    ModuleIO::output_convergence_after_scf(conv_esolver, this->pelec->f_en.etot);

    //! Write Fermi energy
    ModuleIO::output_efermi(conv_esolver, this->pelec->eferm.ef);

    //! Update delta_rho for charge extrapolation
    CE.update_delta_rho(ucell, &(this->chr), &(this->sf));

    //! print out charge density, potential, elf, etc.
	ModuleIO::ctrl_output_fp(ucell, this->pelec, this->pw_big, this->pw_rhod, 
			this->chr, this->solvent, this->Pgrid, istep); 

}

void ESolver_FP::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_FP", "before_scf");

    // if the cell has changed
    if (ucell.cell_parameter_updated)
    {
        // only G-vector and K-vector are changed due to the change of lattice
        // vector FFT grids do not change!!
        this->pw_rho->initgrids(ucell.lat0, ucell.latvec, pw_rho->nx, pw_rho->ny, pw_rho->nz);
        this->pw_rho->collect_local_pw();
        this->pw_rho->collect_uniqgg();

        // if double grid used in USPP, update related quantities in dense grid
        if (PARAM.globalv.double_grid)
        {
            this->pw_rhod->initgrids(ucell.lat0, ucell.latvec, pw_rhod->nx, pw_rhod->ny, pw_rhod->nz);
            this->pw_rhod->collect_local_pw();
            this->pw_rhod->collect_uniqgg();
        }

        // reset local pseudopotentials
        this->locpp.init_vloc(ucell, this->pw_rhod);
        this->locpp.print_vloc(ucell, this->pw_rhod,
                               this->inp_->out_element_info, PARAM.globalv.global_out_dir);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

        // perform symmetry analysis
        if (ModuleSymmetry::Symmetry::symm_flag == 1)
        {
            const int cal_symm_repr[2] = {this->inp_->cal_symm_repr[0], this->inp_->cal_symm_repr[1]};
            ucell.symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running,
                                 this->inp_->symmetry_prec, this->inp_->nspin, this->inp_->calculation, cal_symm_repr);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
        }

        // reset k-points
        kv.set_after_vc(this->inp_->nspin, ucell.G, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT K-POINTS");
    }

    // charge extrapolation
    if (ucell.ionic_position_updated)
    {
        this->CE.update_all_dis(ucell);
        this->CE.extrapolate_charge(&this->Pgrid, ucell, &this->chr, &this->sf,
                                    GlobalV::ofs_running, GlobalV::ofs_warning);
    }

    //! Evaluate the vdW correction once for this ionic configuration.
    this->vdw_result_.reset();
    auto vdw_solver = vdw::make_vdw(ucell, *this->inp_, &(GlobalV::ofs_running));
    if (vdw_solver != nullptr)
    {
        const vdw::VdwRequest request(this->inp_->cal_force, this->inp_->cal_stress);
        this->vdw_result_.reset(new vdw::VdwResult(vdw_solver->evaluate(request)));
        this->pelec->f_en.evdw = this->vdw_result_->energy;
    }
    else
    {
        this->pelec->f_en.evdw = 0.0;
    }

    //! calculate ewald energy
    if (!this->inp_->test_skip_ewald)
    {
        this->pelec->f_en.ewald_energy = H_Ewald_pw::compute_ewald(ucell, this->pw_rhod, this->sf.strucFac, this->inp_->test_energy, GlobalV::ofs_running);
    }

    //! set direction of magnetism, used in non-collinear case 
    unitcell::cal_ux(ucell, this->inp_->nspin);

    //! output the initial charge density
    ModuleIO::write_chg_init(ucell, this->Pgrid, this->chr, this->pelec->eferm, istep,
                             PARAM.globalv.global_out_dir, *this->inp_, PARAM.globalv.two_fermi);

    return;
}

void ESolver_FP::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    //! output charge density in G-space, or if available, kinetic energy density in G-space
    if (this->inp_->out_chg[0] != -1)
    {
        if (iter % this->inp_->out_freq_elec == 0 || iter == this->inp_->scf_nmax || conv_esolver)
        {
            for (int is = 0; is < this->inp_->nspin; is++)
            {
                this->pw_rhod->real2recip(this->chr.rho_save[is], this->chr.rhog_save[is]);
            }
            ModuleIO::write_rhog(PARAM.globalv.global_out_dir + this->inp_->suffix + "-CHARGE-DENSITY.restart",
                                 PARAM.globalv.gamma_only_pw,
                                 this->pw_rhod,
                                 this->inp_->nspin,
                                 ucell.GT,
                                 this->chr.rhog_save,
                                 GlobalV::MY_POOL,
                                 GlobalV::RANK_IN_POOL,
                                 GlobalV::NPROC_IN_POOL);

            if (XC_Functional::get_ked_flag())
            {
                std::vector<std::complex<double>> kin_g_space(this->inp_->nspin * this->chr.ngmc, {0.0, 0.0});
                std::vector<std::complex<double>*> kin_g;
                for (int is = 0; is < this->inp_->nspin; is++)
                {
                    kin_g.push_back(kin_g_space.data() + is * this->chr.ngmc);
                    this->pw_rhod->real2recip(this->chr.kin_r_save[is], kin_g[is]);
                }
                ModuleIO::write_rhog(PARAM.globalv.global_out_dir + this->inp_->suffix + "-TAU-DENSITY.restart",
                                     PARAM.globalv.gamma_only_pw,
                                     this->pw_rhod,
                                     this->inp_->nspin,
                                     ucell.GT,
                                     kin_g.data(),
                                     GlobalV::MY_POOL,
                                     GlobalV::RANK_IN_POOL,
                                     GlobalV::NPROC_IN_POOL);
            }
        }
    }
}

void ESolver_FP::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unit_cell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    // print out the final total energy
    GlobalV::ofs_running << "\n --------------------------------------------" << std::endl;
    GlobalV::ofs_running << std::setprecision(16);
    GlobalV::ofs_running << " !FINAL_ETOT_IS " << this->pelec->f_en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;

}

} // namespace ModuleESolver
