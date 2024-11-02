#include "esolver_ks_pw.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_io/cube_io.h"
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
void ESolver_KS_PW<T, Device>::nscf() {
    ModuleBase::TITLE("ESolver_KS_PW", "nscf");
    ModuleBase::timer::tick("ESolver_KS_PW", "nscf");

    //! 1) before_scf
    const int istep_tmp = 0;
    this->before_scf(istep_tmp);

    //! 2) setup the parameters for diagonalization
    double diag_ethr = PARAM.inp.pw_diag_thr;
    if (diag_ethr - 1e-2 > -1e-5) {
        diag_ethr
            = std::max(1e-13,
                       0.1 * std::min(1e-2, PARAM.inp.scf_thr / PARAM.inp.nelec));
    }
    GlobalV::ofs_running << " PW_DIAG_THR  = " << diag_ethr << std::endl;

    this->hamilt2estates(diag_ethr);

    //! 3) calculate weights/Fermi energies
    this->pelec->calculate_weights();

    GlobalV::ofs_running << "\n End of Band Structure Calculation \n"
                         << std::endl;

    //! 4) print out band energies and weights
    std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "writing band energies");
    const int nspin = PARAM.inp.nspin;
    const int nbands = PARAM.inp.nbands;
    for (int ik = 0; ik < this->kv.get_nks(); ik++) {
        if (nspin == 2) {
            if (ik == 0) {
                GlobalV::ofs_running << " spin up :" << std::endl;
            }
            if (ik == (this->kv.get_nks() / 2)) {
                GlobalV::ofs_running << " spin down :" << std::endl;
            }
        }

        GlobalV::ofs_running
            << " k-points" << ik + 1 << "(" << this->kv.get_nkstot()
            << "): " << this->kv.kvec_c[ik].x << " " << this->kv.kvec_c[ik].y
            << " " << this->kv.kvec_c[ik].z << std::endl;

        for (int ib = 0; ib < nbands; ib++) {
            GlobalV::ofs_running
                << " spin" << this->kv.isk[ik] + 1 << "_final_band " << ib + 1
                << " " << this->pelec->ekb(ik, ib) * ModuleBase::Ry_to_eV << " "
                << this->pelec->wg(ik, ib) * this->kv.get_nks() << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
    std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "writing band energies");

    //! 5) print out band gaps
    if (PARAM.inp.out_bandgap) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "writing band gaps");
        if (!PARAM.globalv.two_fermi) {
            this->pelec->cal_bandgap();
            GlobalV::ofs_running << " E_bandgap "
                                 << this->pelec->bandgap * ModuleBase::Ry_to_eV
                                 << " eV" << std::endl;
        } else {
            this->pelec->cal_bandgap_updw();
            GlobalV::ofs_running
                << " E_bandgap_up "
                << this->pelec->bandgap_up * ModuleBase::Ry_to_eV << " eV"
                << std::endl;
            GlobalV::ofs_running
                << " E_bandgap_dw "
                << this->pelec->bandgap_dw * ModuleBase::Ry_to_eV << " eV"
                << std::endl;
        }
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "writing band gaps");
    }

    //! 6) calculate Wannier functions
    if (PARAM.inp.towannier90) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Wannier functions calculation");
        toWannier90_PW wan(PARAM.inp.out_wannier_mmn,
                           PARAM.inp.out_wannier_amn,
                           PARAM.inp.out_wannier_unk,
                           PARAM.inp.out_wannier_eig,
                           PARAM.inp.out_wannier_wvfn_formatted,
                           PARAM.inp.nnkpfile,
                           PARAM.inp.wannier_spin);

        wan.calculate(this->pelec->ekb,
                      this->pw_wfc,
                      this->pw_big,
                      this->kv,
                      this->psi);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Wannier functions calculation");
    }

    //! 7) calculate Berry phase polarization
    if (berryphase::berry_phase_flag
        && ModuleSymmetry::Symmetry::symm_flag != 1) {
        std::cout << FmtCore::format("\n * * * * * *\n << Start %s.\n", "Berry phase polarization");
        berryphase bp;
        bp.Macroscopic_polarization(this->pw_wfc->npwk_max,
                                    this->psi,
                                    this->pw_rho,
                                    this->pw_wfc,
                                    this->kv);
        std::cout << FmtCore::format(" >> Finish %s.\n * * * * * *\n", "Berry phase polarization");
    }

    /// write potential
    if (PARAM.inp.out_pot == 1 || PARAM.inp.out_pot == 3)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            std::string fn = PARAM.globalv.global_out_dir + "/SPIN" + std::to_string(is + 1) + "_POT.cube";

            ModuleIO::write_vdata_palgrid(GlobalC::Pgrid,
                this->pelec->pot->get_effective_v(is),
                is,
                PARAM.inp.nspin,
                0,
                fn,
                0.0, // efermi
                &(GlobalC::ucell),
                3,  // precision
                0); // out_fermi
        }
    }
    else if (PARAM.inp.out_pot == 2)
    {
        std::string fn = PARAM.globalv.global_out_dir + "/ElecStaticPot.cube";
        ModuleIO::write_elecstat_pot(
#ifdef __MPI
            this->pw_big->bz,
            this->pw_big->nbz,
#endif
            fn,
            0, // istep
            this->pw_rhod,
            this->pelec->charge,
            &(GlobalC::ucell),
            this->pelec->pot->get_fixed_v());
    }

    ModuleBase::timer::tick("ESolver_KS_PW", "nscf");
    return;
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver
