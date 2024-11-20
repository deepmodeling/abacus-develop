#include "esolver_ks_pw.h"

#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/input_conv.h"
#include "module_io/nscf_band.h"
#include "module_io/output_log.h"
#include "module_io/write_dos_pw.h"
#include "module_io/write_istate_info.h"
#include "module_io/write_wfc_pw.h"

#include <iostream>

//--------------temporary----------------------------
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_general/module_ewald/H_Ewald_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"
//-----force-------------------
#include "module_hamilt_pw/hamilt_pwdft/forces.h"
//-----stress------------------
#include "module_hamilt_pw/hamilt_pwdft/stress_pw.h"
//---------------------------------------------------
#include "module_base/memory.h"
#include "module_base/module_device/device.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_io/berryphase.h"
#include "module_io/numerical_basis.h"
#include "module_io/numerical_descriptor.h"
#include "module_io/to_wannier90_pw.h"
#include "module_io/winput.h"
#include "module_io/write_elecstat_pot.h"
#include "module_io/write_wfc_r.h"
#include "module_parameter/parameter.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif
#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include "module_base/formatter.h"

namespace ModuleESolver {

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::init_after_vc(const Input_para& inp, UnitCell& ucell) {
    ModuleBase::TITLE("ESolver_KS_PW", "init_after_vc");
    ModuleBase::timer::tick("ESolver_KS_PW", "init_after_vc");

    ESolver_KS<T, Device>::init_after_vc(inp, ucell);

    if (inp.mdp.md_prec_level == 2) {
        this->pw_wfc->initgrids(ucell.lat0,
                                ucell.latvec,
                                this->pw_rho->nx,
                                this->pw_rho->ny,
                                this->pw_rho->nz);

        this->pw_wfc->initparameters(false,
                                     inp.ecutwfc,
                                     this->kv.get_nks(),
                                     this->kv.kvec_d.data());

#ifdef __MPI
        if (PARAM.inp.pw_seed > 0) {
            MPI_Allreduce(MPI_IN_PLACE,
                          &this->pw_wfc->ggecut,
                          1,
                          MPI_DOUBLE,
                          MPI_MAX,
                          MPI_COMM_WORLD);
        }
        // qianrui add 2021-8-13 to make different kpar parameters can get the
        // same results
#endif

        this->pw_wfc->setuptransform();

        for (int ik = 0; ik < this->kv.get_nks(); ++ik) {
            this->kv.ngk[ik] = this->pw_wfc->npwk[ik];
        }

        this->pw_wfc->collect_local_pw(inp.erf_ecut,
                                       inp.erf_height,
                                       inp.erf_sigma);
        this->already_initpsi = false;

        delete this->pelec;
        this->pelec
            = new elecstate::ElecStatePW<T, Device>(this->pw_wfc,
                                                    &(this->chr),
                                                    (K_Vectors*)(&(this->kv)),
                                                    &ucell,
                                                    &GlobalC::ppcell,
                                                    this->pw_rhod,
                                                    this->pw_rho,
                                                    this->pw_big);

        this->pelec->charge->allocate(PARAM.inp.nspin);

        //! setup cell volume
        this->pelec->omega = ucell.omega;

        delete this->pelec->pot;

        this->pelec->pot = new elecstate::Potential(this->pw_rhod,
                                                    this->pw_rho,
                                                    &ucell,
                                                    &GlobalC::ppcell.vloc,
                                                    &(this->sf),
                                                    &(this->pelec->f_en.etxc),
                                                    &(this->pelec->f_en.vtxc));

        // temporary
        this->Init_GlobalC(inp, ucell, GlobalC::ppcell);
    } else {
        GlobalC::ppcell.init_vnl(ucell, this->pw_rhod);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,
                                     "NON-LOCAL POTENTIAL");

        this->pw_wfc->initgrids(ucell.lat0,
                                ucell.latvec,
                                this->pw_wfc->nx,
                                this->pw_wfc->ny,
                                this->pw_wfc->nz);

        this->pw_wfc->initparameters(false,
                                     PARAM.inp.ecutwfc,
                                     this->kv.get_nks(),
                                     this->kv.kvec_d.data());

        this->pw_wfc->collect_local_pw(inp.erf_ecut,
                                       inp.erf_height,
                                       inp.erf_sigma);

        this->p_wf_init->make_table(this->kv.get_nks(), &this->sf);
    }

#ifdef USE_PAW
    if (PARAM.inp.use_paw) {
        GlobalC::paw_cell.set_libpaw_ecut(PARAM.inp.ecutwfc / 2.0,
                                          PARAM.inp.ecutwfc / 2.0); // in Hartree
        GlobalC::paw_cell.set_libpaw_fft(this->pw_wfc->nx,
                                         this->pw_wfc->ny,
                                         this->pw_wfc->nz,
                                         this->pw_wfc->nx,
                                         this->pw_wfc->ny,
                                         this->pw_wfc->nz,
                                         this->pw_wfc->startz,
                                         this->pw_wfc->numz);

#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0) {
            GlobalC::paw_cell.prepare_paw();
        }
#else
        GlobalC::paw_cell.prepare_paw();
#endif
        GlobalC::paw_cell.set_sij();

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;
#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0) {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < ucell.nat; iat++) {
                GlobalC::paw_cell.set_rhoij(iat,
                                            nrhoijsel[iat],
                                            rhoijselect[iat].size(),
                                            rhoijselect[iat].data(),
                                            rhoijp[iat].data());
            }
        }
#else
        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for (int iat = 0; iat < ucell.nat; iat++) {
            GlobalC::paw_cell.set_rhoij(iat,
                                        nrhoijsel[iat],
                                        rhoijselect[iat].size(),
                                        rhoijselect[iat].data(),
                                        rhoijp[iat].data());
        }
#endif
    }
#endif

    ModuleBase::timer::tick("ESolver_KS_PW", "init_after_vc");
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
