#include "esolver_ks_lcao_tddft.h"

#include "module_io/cal_r_overlap_R.h"
#include "module_io/dipole_io.h"
#include "module_io/td_current_io.h"
#include "module_io/write_HS.h"
#include "module_io/write_HS_R.h"
#include "module_elecstate/elecstate_tools.h"

//--------------temporary----------------------------
#include "source_base/blas_connector.h"
#include "source_base/global_function.h"
#include "source_base/lapack_connector.h"
#include "source_base/scalapack_connector.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_elecstate/module_dm/cal_edm_tddft.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/cal_ux.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "source_hsolver/hsolver_lcao.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//---------------------------------------------------

namespace ModuleESolver
{

template <typename Device>
ESolver_KS_LCAO_TDDFT<Device>::ESolver_KS_LCAO_TDDFT()
{
    classname = "ESolver_rtTDDFT";
    basisname = "LCAO";

    // If the device is GPU, we must open use_tensor and use_lapack
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        use_tensor = true;
        use_lapack = true;
    }
}

template <typename Device>
ESolver_KS_LCAO_TDDFT<Device>::~ESolver_KS_LCAO_TDDFT()
{
    delete psi_laststep;
    if (Hk_laststep != nullptr)
    {
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            delete[] Hk_laststep[ik];
        }
        delete[] Hk_laststep;
    }
    if (Sk_laststep != nullptr)
    {
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            delete[] Sk_laststep[ik];
        }
        delete[] Sk_laststep;
    }
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    // 1) run before_all_runners in ESolver_KS_LCAO
    ESolver_KS_LCAO<std::complex<double>, double>::before_all_runners(ucell, inp);

    // this line should be optimized
    // this->pelec = dynamic_cast<elecstate::ElecStateLCAO_TDDFT*>(this->pelec);
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::hamilt2rho_single(UnitCell& ucell,
                                                          const int istep,
                                                          const int iter,
                                                          const double ethr)
{
    if (PARAM.inp.init_wfc == "file")
    {
        if (istep >= 1)
        {
            module_tddft::Evolve_elec<Device>::solve_psi(istep,
                                                         PARAM.inp.nbands,
                                                         PARAM.globalv.nlocal,
                                                         kv.get_nks(),
                                                         this->p_hamilt,
                                                         this->pv,
                                                         this->psi,
                                                         this->psi_laststep,
                                                         this->Hk_laststep,
                                                         this->Sk_laststep,
                                                         this->pelec->ekb,
                                                         GlobalV::ofs_running,
                                                         td_htype,
                                                         PARAM.inp.propagator,
                                                         use_tensor,
                                                         use_lapack);
            this->weight_dm_rho();
        }
        this->weight_dm_rho();
    }
    else if (istep >= 2)
    {
        module_tddft::Evolve_elec<Device>::solve_psi(istep,
                                                     PARAM.inp.nbands,
                                                     PARAM.globalv.nlocal,
                                                     kv.get_nks(),
                                                     this->p_hamilt,
                                                     this->pv,
                                                     this->psi,
                                                     this->psi_laststep,
                                                     this->Hk_laststep,
                                                     this->Sk_laststep,
                                                     this->pelec->ekb,
                                                     GlobalV::ofs_running,
                                                     td_htype,
                                                     PARAM.inp.propagator,
                                                     use_tensor,
                                                     use_lapack);
        this->weight_dm_rho();
    }
    else
    {
        // reset energy
        this->pelec->f_en.eband = 0.0;
        this->pelec->f_en.demet = 0.0;
        if (this->psi != nullptr)
        {
            bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;
            hsolver::HSolverLCAO<std::complex<double>> hsolver_lcao_obj(&this->pv, PARAM.inp.ks_solver);
            hsolver_lcao_obj.solve(this->p_hamilt, this->psi[0], this->pelec, skip_charge);
        }
    }

    // symmetrize the charge density only for ground state
    if (istep <= 1)
    {
        Symmetry_rho srho;
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            srho.begin(is, this->chr, pw_rho, ucell.symm);
        }
    }

    // (7) calculate delta energy
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::iter_finish(
		UnitCell& ucell, 
		const int istep, 
		int& iter,
		bool& conv_esolver)
{
    // print occupation of each band
    if (iter == 1 && istep <= 2)
    {
        GlobalV::ofs_running << " ---------------------------------------------------------"
                             << std::endl;
        GlobalV::ofs_running << " occupations of electrons" << std::endl;
        GlobalV::ofs_running << " k-point  state   occupation" << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);
        GlobalV::ofs_running << std::left;
        std::setprecision(6);
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                GlobalV::ofs_running << " " << std::setw(9) 
                 << ik+1 << std::setw(8) << ib + 1 
                 << std::setw(12) << this->pelec->wg(ik, ib) << std::endl;
            }
        }
        GlobalV::ofs_running << " ---------------------------------------------------------"
                             << std::endl;
    }

    ESolver_KS_LCAO<std::complex<double>, double>::iter_finish(ucell, istep, iter, conv_esolver);
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::update_pot(UnitCell& ucell, 
		const int istep, 
		const int iter, 
		const bool conv_esolver)
{
    // Calculate new potential according to new Charge Density
    if (!conv_esolver)
    {
        elecstate::cal_ux(ucell);
        this->pelec->pot->update_from_charge(&this->chr, &ucell);
        this->pelec->f_en.descf = this->pelec->cal_delta_escf();
    }
    else
    {
        this->pelec->cal_converged();
    }

    const int nloc = this->pv.nloc;
    const int ncol_nbands = this->pv.ncol_bands;
    const int nrow = this->pv.nrow;
    const int nbands = PARAM.inp.nbands;
    const int nlocal = PARAM.globalv.nlocal;

    // store wfc and Hk laststep
    if (istep >= (PARAM.inp.init_wfc == "file" ? 0 : 1) && conv_esolver)
    {
        if (this->psi_laststep == nullptr)
        {
            int ncol_tmp = 0;
            int nrow_tmp = 0;
#ifdef __MPI
            ncol_tmp = ncol_nbands;
            nrow_tmp = nrow;
#else
            ncol_tmp = nbands;
            nrow_tmp = nlocal;
#endif
            this->psi_laststep = new psi::Psi<std::complex<double>>(kv.get_nks(), ncol_tmp, nrow_tmp, kv.ngk, true);

        }

        // allocate memory for Hk_laststep and Sk_laststep
        if (td_htype == 1)
        {
            // Length of Hk_laststep and Sk_laststep, nlocal * nlocal for global, nloc for local
            const int len_HS = use_tensor && use_lapack ? nlocal * nlocal : nloc;

            if (this->Hk_laststep == nullptr)
            {
                this->Hk_laststep = new std::complex<double>*[kv.get_nks()];
                for (int ik = 0; ik < kv.get_nks(); ++ik)
                {
                    // Allocate memory for Hk_laststep, if (use_tensor && use_lapack), should be global
                    this->Hk_laststep[ik] = new std::complex<double>[len_HS];
                    ModuleBase::GlobalFunc::ZEROS(Hk_laststep[ik], len_HS);
                }
            }
            if (this->Sk_laststep == nullptr)
            {
                this->Sk_laststep = new std::complex<double>*[kv.get_nks()];
                for (int ik = 0; ik < kv.get_nks(); ++ik)
                {
                    // Allocate memory for Sk_laststep, if (use_tensor && use_lapack), should be global
                    this->Sk_laststep[ik] = new std::complex<double>[len_HS];
                    ModuleBase::GlobalFunc::ZEROS(Sk_laststep[ik], len_HS);
                }
            }
        }

        // put information to Hk_laststep and Sk_laststep
        for (int ik = 0; ik < kv.get_nks(); ++ik)
        {
            this->psi->fix_k(ik);
            this->psi_laststep->fix_k(ik);

            // copy the data from psi to psi_laststep
            const int size0 = psi->get_nbands() * psi->get_nbasis();
            for (int index = 0; index < size0; ++index)
            {
                psi_laststep[0].get_pointer()[index] = psi[0].get_pointer()[index];
            }

            // store Hamiltonian
            if (td_htype == 1)
            {
                this->p_hamilt->updateHk(ik);
                hamilt::MatrixBlock<complex<double>> h_mat;
                hamilt::MatrixBlock<complex<double>> s_mat;
                this->p_hamilt->matrix(h_mat, s_mat);

                if (use_tensor && use_lapack)
                {
                    // Gather H and S matrices to root process
#ifdef __MPI
                    int myid = 0;
                    int num_procs = 1;
                    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
                    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

                    Matrix_g<std::complex<double>> h_mat_g; // Global matrix structure
                    Matrix_g<std::complex<double>> s_mat_g; // Global matrix structure

                    // Collect H matrix
                    gatherMatrix(myid, 0, h_mat, h_mat_g);
                    BlasConnector::copy(nlocal * nlocal, h_mat_g.p.get(), 1, Hk_laststep[ik], 1);

                    // Collect S matrix
                    gatherMatrix(myid, 0, s_mat, s_mat_g);
                    BlasConnector::copy(nlocal * nlocal, s_mat_g.p.get(), 1, Sk_laststep[ik], 1);
#endif
                }
                else
                {
                    BlasConnector::copy(nloc, h_mat.p, 1, Hk_laststep[ik], 1);
                    BlasConnector::copy(nloc, s_mat.p, 1, Sk_laststep[ik], 1);
                }
            }
        }

        // calculate energy density matrix for tddft
        if (istep >= (PARAM.inp.init_wfc == "file" ? 0 : 2) && PARAM.inp.td_edm == 0)
        {
            elecstate::cal_edm_tddft(this->pv, this->pelec, this->kv, this->p_hamilt);
        }
    }

    // print "eigen value" for tddft
// it seems uncessary to print out E_ii because the band energies are printed
/*
    if (conv_esolver)
    {
        GlobalV::ofs_running << "----------------------------------------------------------"
                             << std::endl;
        GlobalV::ofs_running << " Print E=<psi_i|H|psi_i> " << std::endl;
        GlobalV::ofs_running << " k-point  state    energy (eV)" << std::endl;
        GlobalV::ofs_running << "----------------------------------------------------------"
                             << std::endl;
        GlobalV::ofs_running << std::setprecision(6);
        GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                GlobalV::ofs_running << " " << std::setw(7) << ik + 1 
                                     << std::setw(7) << ib + 1 
                                     << std::setw(10) << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV 
                                     << std::endl;
            }
        }
    }
*/
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_LCAO_TDDFT", "after_scf");
    ModuleBase::timer::tick("ESolver_LCAO_TDDFT", "after_scf");

    ESolver_KS_LCAO<std::complex<double>, double>::after_scf(ucell, istep, conv_esolver);

    // (1) write dipole information
    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        if (PARAM.inp.out_dipole == 1)
        {
            std::stringstream ss_dipole;
            ss_dipole << PARAM.globalv.global_out_dir << "SPIN" << is + 1 << "_DIPOLE";
            ModuleIO::write_dipole(ucell,
                                   this->chr.rho_save[is],
                                   this->chr.rhopw,
                                   is,
                                   istep,
                                   ss_dipole.str());
        }
    }

     // (2) write current information
    if (TD_Velocity::out_current == true)
    {
        elecstate::DensityMatrix<std::complex<double>, double>* tmp_DM
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM();

        ModuleIO::write_current(ucell,
                                this->gd,
                                istep,
                                this->psi,
                                pelec,
                                kv,
                                two_center_bundle_.overlap_orb.get(),
                                tmp_DM->get_paraV_pointer(),
                                orb_,
                                this->RA);
    }


    ModuleBase::timer::tick("ESolver_LCAO_TDDFT", "after_scf");
}

template <typename Device>
void ESolver_KS_LCAO_TDDFT<Device>::weight_dm_rho()
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

    // calculate Eband energy
    elecstate::calEBand(this->pelec->ekb,this->pelec->wg,this->pelec->f_en);

    // calculate the density matrix
    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");

    auto _pes = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec);
    elecstate::cal_dm_psi(_pes->DM->get_paraV_pointer(), _pes->wg, this->psi[0], *(_pes->DM));
    _pes->DM->cal_DMR();

    // get the real-space charge density
    this->pelec->psiToRho(this->psi[0]);
}

template class ESolver_KS_LCAO_TDDFT<base_device::DEVICE_CPU>;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template class ESolver_KS_LCAO_TDDFT<base_device::DEVICE_GPU>;
#endif

} // namespace ModuleESolver
