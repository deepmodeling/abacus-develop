#include "of_stress_pw.h"

#include "module_base/timer.h"
#include "module_hamilt_general/module_vdw/vdw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/output_log.h"

// Since the kinetic stress of OFDFT is calculated by kinetic functionals in esolver_of.cpp, here we regard it as an
// input variable.
void OF_Stress_PW::cal_stress(ModuleBase::matrix& sigmatot,
                              ModuleBase::matrix& kinetic_stress,
                              UnitCell& ucell,
                              ModuleSymmetry::Symmetry* p_symm,
                              Structure_Factor* p_sf,
                              K_Vectors* p_kv)
{
    ModuleBase::TITLE("OF_Stress_PW", "cal_stress");
    ModuleBase::timer::tick("OF_Stress_PW", "cal_stress");

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

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sigmatot(i, j) = 0.0;
            sigmaxc(i, j) = 0.0;
            sigmahar(i, j) = 0.0;
            // kinetic contribution
            sigmakin(i, j) = kinetic_stress(i, j);
            sigmaloc(i, j) = 0.0;
            sigmanl(i, j) = 0.0;
            sigmaewa(i, j) = 0.0;
            sigmaxcc(i, j) = 0.0;
            sigmavdw(i, j) = 0.0;
        }
    }

    // hartree contribution
    stress_har(sigmahar, this->rhopw, 1, pelec->charge);

    // ewald contribution
    stress_ewa(sigmaewa, this->rhopw, 1);

    // xc contribution: add gradient corrections(non diagonal)
    for (int i = 0; i < 3; i++)
    {
        sigmaxc(i, i) = -(pelec->f_en.etxc - pelec->f_en.vtxc) / ucell.omega;
    }
    stress_gga(sigmaxc, this->rhopw, pelec->charge);

    // local contribution
    stress_loc(sigmaloc, this->rhopw, p_sf, 1, pelec->charge);

    // nlcc
    stress_cc(sigmaxcc, this->rhopw, p_sf, 1, pelec->charge);

    // vdw term
    stress_vdw(sigmavdw, ucell);

    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            sigmatot(ipol, jpol) = sigmakin(ipol, jpol) + sigmahar(ipol, jpol) + sigmanl(ipol, jpol)
                                   + sigmaxc(ipol, jpol) + sigmaxcc(ipol, jpol) + sigmaewa(ipol, jpol)
                                   + sigmaloc(ipol, jpol) + sigmavdw(ipol, jpol);
        }
    }

    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigmatot, ucell);
    }

    bool ry = false;
    ModuleIO::print_stress("TOTAL-STRESS", sigmatot, true, ry);

    if (GlobalV::TEST_STRESS)
    {
        GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
        ModuleIO::print_stress("KINETIC    STRESS", sigmakin, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("LOCAL    STRESS", sigmaloc, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("HARTREE    STRESS", sigmahar, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("NON-LOCAL    STRESS", sigmanl, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("XC    STRESS", sigmaxc, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("EWALD    STRESS", sigmaewa, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("NLCC    STRESS", sigmaxcc, GlobalV::TEST_STRESS, ry);
        ModuleIO::print_stress("TOTAL    STRESS", sigmatot, GlobalV::TEST_STRESS, ry);
    }
    ModuleBase::timer::tick("OF_Stress_PW", "cal_stress");
    return;
}

void OF_Stress_PW::stress_vdw(ModuleBase::matrix& sigma, UnitCell& ucell)
{
    auto vdw_solver = vdw::make_vdw(ucell, INPUT);
    if (vdw_solver != nullptr)
    {
        sigma = vdw_solver->get_stress().to_matrix();
    }
    return;
}
