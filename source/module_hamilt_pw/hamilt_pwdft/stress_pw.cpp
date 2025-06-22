#include "stress_pw.h"

#include "source_base/timer.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::cal_stress(ModuleBase::matrix& sigmatot,
                                           UnitCell& ucell,
                                           const pseudopot_cell_vl& locpp,
                                           const pseudopot_cell_vnl& nlpp,
                                           ModulePW::PW_Basis* rho_basis,
                                           ModuleSymmetry::Symmetry* p_symm,
                                           Structure_Factor* p_sf,
                                           K_Vectors* p_kv,
                                           ModulePW::PW_Basis_K* wfc_basis,
                                           const psi::Psi<complex<FPTYPE>, Device>* d_psi_in)
{
    ModuleBase::TITLE("Stress_PW", "cal_stress");
    ModuleBase::timer::tick("Stress_PW", "cal_stress");

    // total stress
    sigmatot.create(3, 3);
    ModuleBase::matrix sigmaxc;
    // exchange-correlation stress
    sigmaxc.create(3, 3);
    // hartree stress
    ModuleBase::matrix sigmahar;
    sigmahar.create(3, 3);
    // electron kinetic stress
    ModuleBase::matrix sigmakin;
    sigmakin.create(3, 3);
    // local pseudopotential stress
    ModuleBase::matrix sigmaloc;
    sigmaloc.create(3, 3);
    // non-local pseudopotential stress
    ModuleBase::matrix sigmanl;
    sigmanl.create(3, 3);
    // Ewald stress
    ModuleBase::matrix sigmaewa;
    sigmaewa.create(3, 3);
    // non-linear core correction stress
    ModuleBase::matrix sigmaxcc;
    sigmaxcc.create(3, 3);
    // vdw stress
    ModuleBase::matrix sigmavdw;
    sigmavdw.create(3, 3);
    // DFT+U and DeltaSpin stress
    ModuleBase::matrix sigmaonsite;
    sigmaonsite.create(3, 3);
    // EXX PW stress
    ModuleBase::matrix sigmaexx;
    sigmaexx.create(3, 3);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sigmatot(i, j) = 0.0;
            sigmaxc(i, j) = 0.0;
            sigmahar(i, j) = 0.0;
            sigmakin(i, j) = 0.0;
            sigmaloc(i, j) = 0.0;
            sigmanl(i, j) = 0.0;
            sigmaewa(i, j) = 0.0;
            sigmaxcc(i, j) = 0.0;
            sigmavdw(i, j) = 0.0;
            sigmaonsite(i, j) = 0.0;
            sigmaexx(i, j) = 0.0;
        }
    }

    // kinetic contribution
    this->stress_kin(sigmakin, this->pelec->wg, p_symm, p_kv, wfc_basis, ucell, d_psi_in);

    // hartree contribution
    this->stress_har(ucell, sigmahar, rho_basis, 1, pelec->charge);

    // ewald contribution
    this->stress_ewa(ucell, sigmaewa, rho_basis, 1);

    // xc contribution: add gradient corrections(non diagonal)
    for (int i = 0; i < 3; i++)
    {
        sigmaxc(i, i) = -(pelec->f_en.etxc - pelec->f_en.vtxc) / ucell.omega;
    }
    this->stress_gga(ucell, sigmaxc, rho_basis, pelec->charge);
    if (XC_Functional::get_ked_flag())
    {
        this->stress_mgga(ucell,
                          sigmaxc,
                          this->pelec->wg,
                          this->pelec->pot->get_effective_vofk(),
                          pelec->charge,
                          p_kv,
                          wfc_basis,
                          d_psi_in);
    }

    // local contribution
    this->stress_loc(ucell, sigmaloc, rho_basis, locpp.vloc, p_sf, 1, pelec->charge);

    // nlcc
    this->stress_cc(sigmaxcc, rho_basis, ucell, p_sf, 1, locpp.numeric, pelec->charge);

    // nonlocal
    this->stress_nl(sigmanl, this->pelec->wg, this->pelec->ekb, p_sf, p_kv, p_symm, wfc_basis, d_psi_in, nlpp, ucell);

    // add US term from augmentation charge derivatives
    if (PARAM.globalv.use_uspp)
    {
        this->stress_us(sigmanl, rho_basis, nlpp, ucell);
    }

    // vdw term
    stress_vdw(sigmavdw, ucell);

    // DFT+U and DeltaSpin stress
    if (PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
    {
        this->stress_onsite(sigmaonsite, this->pelec->wg, wfc_basis, ucell, d_psi_in, p_symm);
    }

    // EXX PW stress
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        this->stress_exx(sigmaexx, this->pelec->wg, rho_basis, wfc_basis, p_kv, d_psi_in, ucell);
    }


    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            sigmatot(ipol, jpol) = sigmakin(ipol, jpol) + sigmahar(ipol, jpol) + sigmanl(ipol, jpol)
                                   + sigmaxc(ipol, jpol) + sigmaxcc(ipol, jpol) + sigmaewa(ipol, jpol)
                                   + sigmaloc(ipol, jpol) + sigmavdw(ipol, jpol) + sigmaonsite(ipol, jpol)
                                   + sigmaexx(ipol, jpol);
        }
    }

    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigmatot, ucell.lat);
    }

    const bool ry = false;
    const bool screen = PARAM.inp.test_stress;
    ModuleIO::print_stress("TOTAL-STRESS", sigmatot, screen, ry, GlobalV::ofs_running);

    if (screen)
    {
        GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
        ModuleIO::print_stress("KINETIC    STRESS", sigmakin, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("LOCAL    STRESS", sigmaloc, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("HARTREE    STRESS", sigmahar, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("NON-LOCAL    STRESS", sigmanl, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("XC    STRESS", sigmaxc, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("EWALD    STRESS", sigmaewa, screen, ry, GlobalV::ofs_running);
        ModuleIO::print_stress("NLCC    STRESS", sigmaxcc, screen, ry, GlobalV::ofs_running);
        if (PARAM.inp.dft_plus_u || PARAM.inp.sc_mag_switch)
        {
            ModuleIO::print_stress("ONSITE    STRESS", sigmaonsite, screen, ry, GlobalV::ofs_running);
        }
        if (GlobalC::exx_info.info_global.cal_exx)
        {
            ModuleIO::print_stress("EXX    STRESS", sigmaexx, screen, ry, GlobalV::ofs_running);
        }
        ModuleIO::print_stress("TOTAL    STRESS", sigmatot, screen, ry, GlobalV::ofs_running);
    }
    ModuleBase::timer::tick("Stress_PW", "cal_stress");
    return;
}

template <typename FPTYPE, typename Device>
void Stress_PW<FPTYPE, Device>::stress_vdw(ModuleBase::matrix& sigma, UnitCell& ucell)
{
    auto vdw_solver = vdw::make_vdw(ucell, PARAM.inp);
    if (vdw_solver != nullptr)
    {
        sigma = vdw_solver->get_stress().to_matrix();
    }
    return;
}

template class Stress_PW<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_PW<double, base_device::DEVICE_GPU>;
#endif
