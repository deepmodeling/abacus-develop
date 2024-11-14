#include "sto_stress_pw.h"

#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_io/output_log.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_nonlocal_tools.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_kin_tools.h"

template <typename FPTYPE, typename Device>
void Sto_Stress_PW<FPTYPE, Device>::cal_stress(ModuleBase::matrix& sigmatot,
                                               const elecstate::ElecState& elec,
                                               ModulePW::PW_Basis* rho_basis,
                                               ModuleSymmetry::Symmetry* p_symm,
                                               Structure_Factor* p_sf,
                                               K_Vectors* p_kv,
                                               ModulePW::PW_Basis_K* wfc_basis,
                                               const psi::Psi<complex<FPTYPE>, Device>& psi_in,
                                               const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf,
                                               const Charge* const chr,
                                               pseudopot_cell_vnl* nlpp,
                                               const UnitCell& ucell_in)
{
    ModuleBase::TITLE("Sto_Stress_PW", "cal_stress");
    ModuleBase::timer::tick("Sto_Stress_PW", "cal_stress");
    const ModuleBase::matrix& wg = elec.wg;
    this->ucell = &ucell_in;
    sigmatot.create(3, 3);
    ModuleBase::matrix sigmaxc(3, 3);
    ModuleBase::matrix sigmahar(3, 3);
    ModuleBase::matrix sigmakin(3, 3);
    ModuleBase::matrix sigmaloc(3, 3);
    ModuleBase::matrix sigmanl(3, 3);
    ModuleBase::matrix sigmaewa(3, 3);
    ModuleBase::matrix sigmaxcc(3, 3);

    // kinetic contribution
    this->sto_stress_kin(sigmakin, wg, p_symm, p_kv, wfc_basis, psi_in, stowf);

    // hartree contribution
    this->stress_har(sigmahar, rho_basis, true, chr);

    // ewald contribution
    this->stress_ewa(sigmaewa, rho_basis, true);

    // xc contribution: add gradient corrections(non diagonal)
    for (int i = 0; i < 3; ++i)
    {
        sigmaxc(i, i) = -(elec.f_en.etxc - elec.f_en.vtxc) / this->ucell->omega;
    }
    this->stress_gga(sigmaxc, rho_basis, chr);

    // local contribution
    this->stress_loc(sigmaloc, rho_basis, p_sf, true, chr);

    // nlcc
    this->stress_cc(sigmaxcc, rho_basis, p_sf, true, chr);

    // nonlocal
    this->sto_stress_nl(sigmanl, wg, p_sf, p_symm, p_kv, wfc_basis, *nlpp, ucell_in, psi_in, stowf);

    for (int ipol = 0; ipol < 3; ++ipol)
    {
        for (int jpol = 0; jpol < 3; ++jpol)
        {
            sigmatot(ipol, jpol) = sigmakin(ipol, jpol) + sigmahar(ipol, jpol) + sigmanl(ipol, jpol)
                                   + sigmaxc(ipol, jpol) + sigmaxcc(ipol, jpol) + sigmaewa(ipol, jpol)
                                   + sigmaloc(ipol, jpol);
        }
    }

    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigmatot, this->ucell->lat);
    }

    bool ry = false;
    ModuleIO::print_stress("TOTAL-STRESS", sigmatot, true, ry);

    if (PARAM.inp.test_stress)
    {
        ry = true;
        GlobalV::ofs_running << "\n PARTS OF STRESS: " << std::endl;
        GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
        GlobalV::ofs_running << std::setiosflags(std::ios::fixed) << std::setprecision(10) << std::endl;
        ModuleIO::print_stress("KINETIC    STRESS", sigmakin, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("LOCAL    STRESS", sigmaloc, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("HARTREE    STRESS", sigmahar, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("NON-LOCAL    STRESS", sigmanl, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("XC    STRESS", sigmaxc, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("EWALD    STRESS", sigmaewa, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("NLCC    STRESS", sigmaxcc, PARAM.inp.test_stress, ry);
        ModuleIO::print_stress("TOTAL    STRESS", sigmatot, PARAM.inp.test_stress, ry);
    }
    ModuleBase::timer::tick("Sto_Stress_PW", "cal_stress");
    return;
}

template <typename FPTYPE, typename Device>
void Sto_Stress_PW<FPTYPE, Device>::sto_stress_kin(ModuleBase::matrix& sigma,
                                                   const ModuleBase::matrix& wg,
                                                   ModuleSymmetry::Symmetry* p_symm,
                                                   K_Vectors* p_kv,
                                                   ModulePW::PW_Basis_K* wfc_basis,
                                                   const psi::Psi<std::complex<FPTYPE>, Device>& psi,
                                                   const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf)
{
    ModuleBase::TITLE("Sto_Stress_PW", "stress_kin");
    ModuleBase::timer::tick("Sto_Stress_PW", "stress_kin");

    int nksbands = psi.get_nbands();
    if (GlobalV::MY_STOGROUP != 0)
    {
        nksbands = 0;
    }
    
    hamilt::FS_Kin_tools<FPTYPE, Device> kin_tool(*this->ucell, p_kv, wfc_basis, wg);

    for (int ik = 0; ik < wfc_basis->nks; ++ik)
    {
        const int stobands = stowf.nchip[ik];
        psi.fix_k(ik);
        stowf.shchi->fix_k(ik);
        
        kin_tool.cal_gk(ik);

        kin_tool.cal_stress_kin(ik, nksbands, true, psi.get_pointer());
        kin_tool.cal_stress_kin(ik, stobands, false, stowf.shchi->get_pointer());
    }

    kin_tool.symmetrize_stress(p_symm, sigma);
    ModuleBase::timer::tick("Sto_Stress_PW", "stress_kin");

    return;
}

template <typename FPTYPE, typename Device>
void Sto_Stress_PW<FPTYPE, Device>::sto_stress_nl(ModuleBase::matrix& sigma,
                                                  const ModuleBase::matrix& wg,
                                                  Structure_Factor* p_sf,
                                                  ModuleSymmetry::Symmetry* p_symm,
                                                  K_Vectors* p_kv,
                                                  ModulePW::PW_Basis_K* wfc_basis,
                                                  const pseudopot_cell_vnl& nlpp,
                                                  const UnitCell& ucell,
                                                  const psi::Psi<std::complex<FPTYPE>, Device>& psi_in,
                                                  const Stochastic_WF<std::complex<FPTYPE>, Device>& stowf)
{
    ModuleBase::TITLE("Sto_Stress_Func", "stres_nl");
    const int nkb = nlpp.nkb;
    if (nkb == 0)
    {
        return;
    }

    ModuleBase::timer::tick("Sto_Stress_Func", "stres_nl");

    int* nchip = stowf.nchip;
    const int npwx = wfc_basis->npwk_max;
    int nksbands = psi_in.get_nbands();
    if (GlobalV::MY_STOGROUP != 0)
    {
        nksbands = 0;
    }

    // allocate memory for the stress
    FPTYPE* stress_device = nullptr;
    resmem_var_op()(this->ctx, stress_device, 9);
    setmem_var_op()(this->ctx, stress_device, 0, 9);
    std::vector<FPTYPE> sigmanlc(9, 0.0);

    hamilt::FS_Nonlocal_tools<FPTYPE, Device> nl_tools(&nlpp, &ucell, p_kv, wfc_basis, p_sf, wg, nullptr);

    for (int ik = 0; ik < p_kv->get_nks(); ik++)
    {
        const int nstobands = nchip[ik];
        const int max_nbands = stowf.shchi->get_nbands() + nksbands;
        const int npw = wfc_basis->npwk[ik];
        psi_in.fix_k(ik);
        stowf.shchi->fix_k(ik);
        nl_tools.cal_vkb(ik, max_nbands);
        // calculate becp = <psi|beta> for all beta functions
        nl_tools.cal_becp(ik, nksbands, psi_in.get_pointer(), 0);
        nl_tools.cal_becp(ik, nstobands, stowf.shchi->get_pointer(), nksbands);
        nl_tools.reduce_pool_becp(max_nbands);
        // calculate dbecp = <psi|d(beta)/dR> for all beta functions
        // calculate stress = \sum <psi|d(beta_j)/dR> * <psi|beta_i> * D_{ij}
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int jpol = 0; jpol <= ipol; jpol++)
            {
                nl_tools.cal_vkb_deri_s(ik, max_nbands, ipol, jpol);
                nl_tools.cal_dbecp_s(ik, nksbands, psi_in.get_pointer(), 0);
                nl_tools.cal_dbecp_s(ik, nstobands, stowf.shchi->get_pointer(), nksbands);
                nl_tools.cal_stress(ik, nksbands, true, ipol, jpol, stress_device, 0);
                nl_tools.cal_stress(ik, nstobands, false, ipol, jpol, stress_device, nksbands);
            }
        }

    }

    // transfer stress from device to host
    syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, sigmanlc.data(), stress_device, 9);
    delmem_var_op()(this->ctx, stress_device);
    // sum up forcenl from all processors
    for (int l = 0; l < 3; l++)
    {
        for (int m = 0; m < 3; m++)
        {
            if (m > l)
            {
                sigmanlc[l * 3 + m] = sigmanlc[m * 3 + l];
            }
        }
    }
    // sum up forcenl from all processors
    Parallel_Reduce::reduce_all(sigmanlc.data(), 9);

    for (int ipol = 0; ipol < 3; ++ipol)
    {
        for (int jpol = 0; jpol < 3; ++jpol)
        {
            sigma(ipol, jpol) = sigmanlc[ipol * 3 + jpol] / ucell.omega;
        }
    }
    // do symmetry
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigma, ucell.lat);
    }

    ModuleBase::timer::tick("Sto_Stress_Func", "stres_nl");
    return;
}


template class Sto_Stress_PW<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Sto_Stress_PW<double, base_device::DEVICE_GPU>;
#endif