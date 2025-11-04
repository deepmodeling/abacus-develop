#include "esolver_ks_lcao_tddft.h"

//----------------IO-----------------
#include "source_io/ctrl_output_td.h"
#include "source_io/dipole_io.h"
#include "source_io/output_log.h"
#include "source_io/read_wfc_nao.h"
#include "source_io/td_current_io.h"
//------LCAO HSolver ElecState-------
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_charge/symmetry_rho.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_estate/module_dm/cal_edm_tddft.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_hsolver/hsolver_lcao.h"
#include "source_lcao/module_rt/evolve_elec.h"
#include "source_lcao/rho_tau_lcao.h"

namespace ModuleESolver
{

template <typename TR, typename Device>
ESolver_KS_LCAO_TDDFT<TR, Device>::ESolver_KS_LCAO_TDDFT()
{
    this->classname = "ESolver_rtTDDFT";
    this->basisname = "LCAO";

    // If the device is GPU, we must open use_tensor and use_lapack
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        use_tensor = true;
        use_lapack = true;
    }
}

template <typename TR, typename Device>
ESolver_KS_LCAO_TDDFT<TR, Device>::~ESolver_KS_LCAO_TDDFT()
{
    //*************************************************
    // Do not add any code in this destructor function
    //*************************************************
    delete psi_laststep;

    if (td_p != nullptr)
    {
        delete td_p;
    }
    TD_info::td_vel_op = nullptr;
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    // Run before_all_runners in ESolver_KS_LCAO
    ESolver_KS_LCAO<std::complex<double>, TR>::before_all_runners(ucell, inp);

    td_p = new TD_info(&ucell);
    TD_info::td_vel_op = td_p;
    totstep += TD_info::estep_shift;

    if (PARAM.inp.init_wfc == "file")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
                                    this->pv,
                                    *(this->psi),
                                    this->pelec->ekb,
                                    this->pelec->wg,
                                    this->kv.ik2iktot,
                                    this->kv.get_nkstot(),
                                    PARAM.inp.nspin,
                                    0,
                                    TD_info::estep_shift))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO_TDDFT", "Read electronic wavefunction from file failed!");
        }
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_LCAO_TDDFT", "runner");
    ModuleBase::timer::tick(this->classname, "runner");

    //----------------------------------------------------------------
    // 1) before_scf (electronic iteration loops)
    //----------------------------------------------------------------
    this->before_scf(ucell, istep); // From ESolver_KS_LCAO
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT SCF");

    // Initialize velocity operator for current calculation
    if (PARAM.inp.td_stype != 1 && TD_info::out_current)
    {
        // initialize the velocity operator
        velocity_mat = new Velocity_op<TR>(&ucell,
                                           &(this->gd),
                                           &this->pv,
                                           this->orb_,
                                           this->two_center_bundle_.overlap_orb.get());
        // calculate velocity operator
        velocity_mat->calculate_grad_term();
        velocity_mat->calculate_vcomm_r();
    }
    int estep_max = (istep == 0 && !PARAM.inp.mdp.md_restart) ? 1 : PARAM.inp.estep_per_md;
    if (PARAM.inp.mdp.md_nstep == 0)
    {
        estep_max = PARAM.inp.estep_per_md + 1;
    }
    // int estep_max = PARAM.inp.estep_per_md;
    for (int estep = 0; estep < estep_max; estep++)
    {
        // calculate total time step
        this->totstep++;
        this->print_step();
        // update At
        if (PARAM.inp.td_stype > 0)
        {
            elecstate::H_TDDFT_pw::update_At();
            td_p->cal_cart_At(elecstate::H_TDDFT_pw::At);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Ax(t)", TD_info::cart_At[0]);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Ay(t)", TD_info::cart_At[1]);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Cartesian vector potential Az(t)", TD_info::cart_At[2]);
        }

        if (estep != 0)
        {
            this->CE.update_all_dis(ucell);
            this->CE.extrapolate_charge(&this->Pgrid,
                                        ucell,
                                        &this->chr,
                                        &this->sf,
                                        GlobalV::ofs_running,
                                        GlobalV::ofs_warning);
            // need to test if correct when estep>0
            this->pelec->init_scf(totstep, ucell, this->Pgrid, this->sf.strucFac, this->locpp.numeric, ucell.symm);

            if (totstep <= PARAM.inp.td_tend + 1)
            {
                TD_info::evolve_once = true;
            }
        }
        //----------------------------------------------------------------
        // 2) SCF iterations
        //----------------------------------------------------------------
        bool conv_esolver = false;
        this->niter = this->maxniter;
        this->diag_ethr = PARAM.inp.pw_diag_thr;
        for (int iter = 1; iter <= this->maxniter; ++iter)
        {
            ModuleIO::write_head_td(GlobalV::ofs_running, istep, totstep, iter, this->basisname);

            // 3) Initialization of SCF iterations
            this->iter_init(ucell, totstep, iter); // From ESolver_KS_LCAO

            // 4) Use Hamiltonian to obtain charge density
            this->hamilt2rho(ucell, totstep, iter, this->diag_ethr); // From ESolver_KS

            // 5) Finish SCF iterations
            this->iter_finish(ucell, totstep, estep, estep_max, iter, conv_esolver);

            // 6) Check convergence
            if (conv_esolver || this->oscillate_esolver)
            {
                this->niter = iter;
                if (this->oscillate_esolver)
                {
                    std::cout << " !! Density oscillation is found, STOP HERE !!" << std::endl;
                }
                break;
            }
        } // end SCF iterations

        //----------------------------------------------------------------
        // 7) after_scf
        //----------------------------------------------------------------
        this->after_scf(ucell, totstep, conv_esolver);
        if (!restart_done && PARAM.inp.mdp.md_restart)
        {
            restart_done = true;
            estep += TD_info::estep_shift % PARAM.inp.estep_per_md;
            if (estep == 0)
            {
                break;
            }
            if (PARAM.inp.mdp.md_nstep != 0)
            {
                estep -= 1;
            }
        }
    }

    if (PARAM.inp.td_stype != 1 && TD_info::out_current)
    {
        delete velocity_mat;
    }

    ModuleBase::timer::tick(this->classname, "runner");
    return;
}

// Output electronic step information
template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::print_step()
{
    std::cout << " -------------------------------------------" << std::endl;
    std::cout << " STEP OF ELECTRON EVOLVE : " << unsigned(totstep) << std::endl;
    std::cout << " -------------------------------------------" << std::endl;
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::hamilt2rho_single(UnitCell& ucell,
                                                          const int istep,
                                                          const int iter,
                                                          const double ethr)
{
    if (PARAM.inp.init_wfc == "file")
    {
        if (istep >= TD_info::estep_shift + 1)
        {
            module_rt::Evolve_elec<Device>::solve_psi(istep,
                                                      PARAM.inp.nbands,
                                                      PARAM.globalv.nlocal,
                                                      this->kv.get_nks(),
                                                      this->p_hamilt,
                                                      this->pv,
                                                      this->psi,
                                                      this->psi_laststep,
                                                      this->Hk_laststep,
                                                      this->Sk_laststep,
                                                      this->pelec->ekb,
                                                      GlobalV::ofs_running,
                                                      PARAM.inp.propagator,
                                                      use_tensor,
                                                      use_lapack);
        }
        this->weight_dm_rho(ucell);
    }
    else if (istep >= 1)
    {
        module_rt::Evolve_elec<Device>::solve_psi(istep,
                                                  PARAM.inp.nbands,
                                                  PARAM.globalv.nlocal,
                                                  this->kv.get_nks(),
                                                  this->p_hamilt,
                                                  this->pv,
                                                  this->psi,
                                                  this->psi_laststep,
                                                  this->Hk_laststep,
                                                  this->Sk_laststep,
                                                  this->pelec->ekb,
                                                  GlobalV::ofs_running,
                                                  PARAM.inp.propagator,
                                                  use_tensor,
                                                  use_lapack);
        this->weight_dm_rho(ucell);
    }
    else
    {
        // For the first step, do normal SCF calculation to get initial state
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;
            hsolver::HSolverLCAO<std::complex<double>> hsolver_lcao_obj(&this->pv, PARAM.inp.ks_solver);
            hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, this->chr, PARAM.inp.nspin, skip_charge);
        }
    }

    // Symmetrize the charge density only for ground state
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, this->chr, this->pw_rho, ucell.symm);
        }
    }

    // Calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::iter_finish(UnitCell& ucell,
                                                    const int istep,
                                                    const int estep,
                                                    const int estep_max,
                                                    int& iter,
                                                    bool& conv_esolver)
{
    // Print occupation of each band
    if (iter == 1 && istep <= 2)
    {
        GlobalV::ofs_running << " k-point  State   Occupations" << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);
        GlobalV::ofs_running << std::left;
        std::setprecision(6);
        for (int ik = 0; ik < this->kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                GlobalV::ofs_running << " " << std::setw(9) << ik + 1 << std::setw(8) << ib + 1 << std::setw(12)
                                     << this->pelec->wg(ik, ib) << std::endl;
            }
        }
        GlobalV::ofs_running << std::endl;
    }

    ESolver_KS_LCAO<std::complex<double>, TR>::iter_finish(ucell, istep, iter, conv_esolver);

    // Store wave function, Hamiltonian and Overlap matrix, to be used in next time step
    this->store_h_s_psi(ucell, istep, iter, conv_esolver);

    // Calculate energy-density matrix for RT-TDDFT
    if (conv_esolver && estep == estep_max - 1 && istep >= (PARAM.inp.init_wfc == "file" ? 0 : 1)
        && PARAM.inp.td_edm == 0)
    {
        elecstate::cal_edm_tddft(this->pv, this->pelec, this->kv, this->p_hamilt);
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::store_h_s_psi(UnitCell& ucell,
                                                      const int istep,
                                                      const int iter,
                                                      const bool conv_esolver)
{
    const int nbands = PARAM.inp.nbands;
    const int nlocal = PARAM.globalv.nlocal;

    // Store wave function, Hamiltonian and Overlap matrix
    if (conv_esolver)
    {
        if (this->psi_laststep == nullptr)
        {
            this->psi_laststep = new psi::Psi<std::complex<double>>(this->kv.get_nks(),
#ifdef __MPI
                                                                    this->pv.ncol_bands,
                                                                    this->pv.nrow,
#else
                                                                    nbands,
                                                                    nlocal,
#endif
                                                                    this->kv.ngk,
                                                                    true);
        }

        // Length of Hk_laststep and Sk_laststep, nlocal * nlocal for global, nloc for local
        const int len_HS_ik = use_tensor && use_lapack ? nlocal * nlocal : this->pv.nloc;
        const int len_HS_all = this->kv.get_nks() * len_HS_ik;

        // Allocate memory for Hk_laststep, if (use_tensor && use_lapack), should be global
        if (this->Hk_laststep.NumElements() != len_HS_all)
        {
            this->Hk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                           ct_device_type_hs,
                                           ct::TensorShape({this->kv.get_nks(), len_HS_ik}));
            this->Hk_laststep.zero();
        }

        // Allocate memory for Sk_laststep, if (use_tensor && use_lapack), should be global
        if (this->Sk_laststep.NumElements() != len_HS_all)
        {
            this->Sk_laststep = ct::Tensor(ct::DataType::DT_COMPLEX_DOUBLE,
                                           ct_device_type_hs,
                                           ct::TensorShape({this->kv.get_nks(), len_HS_ik}));
            this->Sk_laststep.zero();
        }

        // Put information into psi_laststep, Hk_laststep and Sk_laststep
        for (int ik = 0; ik < this->kv.get_nks(); ++ik)
        {
            this->psi->fix_k(ik);
            this->psi_laststep->fix_k(ik);

            // Copy data from psi to psi_laststep at k-point ik
            const int len_psi_ik = this->psi->get_nbands() * this->psi->get_nbasis();
            for (int index = 0; index < len_psi_ik; ++index)
            {
                psi_laststep[0].get_pointer()[index] = this->psi[0].get_pointer()[index];
            }

            // Get H and S matrices at k-point ik
            this->p_hamilt->updateHk(ik);
            hamilt::MatrixBlock<std::complex<double>> h_mat;
            hamilt::MatrixBlock<std::complex<double>> s_mat;
            this->p_hamilt->matrix(h_mat, s_mat);

            // Store H and S matrices to Hk_laststep and Sk_laststep
            if (use_tensor && use_lapack)
            {
                // Gather H and S matrices to root process
#ifdef __MPI
                int myid = 0;
                int num_procs = 1;
                MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

                // Global matrix structure
                module_rt::Matrix_g<std::complex<double>> h_mat_g;
                module_rt::Matrix_g<std::complex<double>> s_mat_g;

                // Collect H matrix
                module_rt::gatherMatrix(myid, 0, h_mat, h_mat_g);
                BlasConnector::copy(len_HS_ik,
                                    h_mat_g.p.get(),
                                    1,
                                    this->Hk_laststep.template data<std::complex<double>>() + ik * len_HS_ik,
                                    1);

                // Collect S matrix
                module_rt::gatherMatrix(myid, 0, s_mat, s_mat_g);
                BlasConnector::copy(len_HS_ik,
                                    s_mat_g.p.get(),
                                    1,
                                    this->Sk_laststep.template data<std::complex<double>>() + ik * len_HS_ik,
                                    1);
#endif
            }
            else
            {
                BlasConnector::copy(len_HS_ik,
                                    h_mat.p,
                                    1,
                                    this->Hk_laststep.template data<std::complex<double>>() + ik * len_HS_ik,
                                    1);
                BlasConnector::copy(len_HS_ik,
                                    s_mat.p,
                                    1,
                                    this->Sk_laststep.template data<std::complex<double>>() + ik * len_HS_ik,
                                    1);
            }
        }
    }
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_LCAO_TDDFT", "after_scf");
    ModuleBase::timer::tick(this->classname, "after_scf");

    ESolver_KS_LCAO<std::complex<double>, TR>::after_scf(ucell, istep, conv_esolver);

    // Output energy for sub-loop (electronic step)
    std::cout << " Potential (Ry): " << std::setprecision(15) << this->pelec->f_en.etot << std::endl;

    // Output dipole, current, etc.
    ModuleIO::ctrl_output_td<TR>(ucell,
                                 this->chr.rho_save,
                                 this->chr.rhopw,
                                 istep,
                                 this->psi,
                                 this->pelec,
                                 this->kv,
                                 this->two_center_bundle_.overlap_orb.get(),
                                 &this->pv,
                                 this->orb_,
                                 this->velocity_mat,
                                 this->RA,
                                 this->td_p);

    ModuleBase::timer::tick(this->classname, "after_scf");
}

template <typename TR, typename Device>
void ESolver_KS_LCAO_TDDFT<TR, Device>::weight_dm_rho(const UnitCell& ucell)
{
    if (PARAM.inp.ocp == 1)
    {
        elecstate::fixed_weights(PARAM.inp.ocp_kb,
                                 PARAM.inp.nbands,
                                 PARAM.inp.nelec,
                                 this->pelec->klist,
                                 this->pelec->wg,
                                 this->pelec->skip_weights);
    }

    // Calculate Eband energy
    elecstate::calEBand(this->pelec->ekb, this->pelec->wg, this->pelec->f_en);

    // Calculate the density matrix
    auto _pes = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
    elecstate::cal_dm_psi(_pes->DM->get_paraV_pointer(), _pes->wg, this->psi[0], *(_pes->DM));
    if (PARAM.inp.td_stype == 2)
    {
        _pes->DM->cal_DMR_td(ucell, TD_info::cart_At);
    }
    else
    {
        _pes->DM->cal_DMR();
    }

    // Get the real-space charge density, mohan add 2025-10-24
    LCAO_domain::dm2rho(_pes->DM->get_DMR_vector(), PARAM.inp.nspin, &this->chr);
}

template class ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_CPU>;
template class ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template class ESolver_KS_LCAO_TDDFT<double, base_device::DEVICE_GPU>;
template class ESolver_KS_LCAO_TDDFT<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
