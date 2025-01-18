#include "module_base/module_device/device.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_lcao/module_deltaspin/spin_constrain.h"
#include "stress_func.h"
// calculate the nonlocal pseudopotential stress in PW
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_onsite(ModuleBase::matrix& sigma,
                                            const ModuleBase::matrix& wg,
                                            const ModulePW::PW_Basis_K* wfc_basis,
                                            const UnitCell& ucell_in,
                                            const psi::Psi<complex<FPTYPE>, Device>* psi_in,
                                            ModuleSymmetry::Symmetry* p_symm)
{
    ModuleBase::TITLE("Stress_Func", "stress_onsite");
    if(psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    ModuleBase::timer::tick("Stress_Func", "stress_onsite");

    FPTYPE* stress_device = nullptr;
    resmem_var_op()(stress_device, 9);
    setmem_var_op()(stress_device, 0, 9);
    std::vector<FPTYPE> sigma_onsite(9, 0.0);

    auto* onsite_p = projectors::OnsiteProjector<FPTYPE, Device>::get_instance();

    const int nks = wfc_basis->nks;
    for (int ik = 0; ik < nks; ik++) // loop k points
    {
        // skip zero weights to speed up
        int nbands_occ = wg.nc;
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int npm = nbands_occ;

        // calculate becp = <psi|beta> for all beta functions
        onsite_p->get_fs_tools()->cal_becp(ik, npm);
        // calculate dbecp = <psi|d(beta)/dR> for all beta functions
        // calculate stress = \sum <psi|d(beta_j)/dR> * <psi|beta_i> * D_{ij}
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int jpol = 0; jpol <= ipol; jpol++)
            {
                FPTYPE* stress_device_tmp = stress_device + (ipol * 3 + jpol);
                onsite_p->get_fs_tools()->cal_dbecp_s(ik, npm, ipol, jpol);
                if(PARAM.inp.dft_plus_u)
                {
                    auto* dftu = ModuleDFTU::DFTU::get_instance();
                    onsite_p->get_fs_tools()->cal_stress_dftu(ik, npm, stress_device_tmp, dftu->orbital_corr.data(), dftu->get_eff_pot_pw(0), dftu->get_size_eff_pot_pw(), wg.c);
                }
                if(PARAM.inp.sc_mag_switch)
                {
                    spinconstrain::SpinConstrain<std::complex<double>>& sc = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
                    const std::vector<ModuleBase::Vector3<double>>& lambda = sc.get_sc_lambda();
                    onsite_p->get_fs_tools()->cal_stress_dspin(ik, npm, stress_device_tmp, lambda.data(), wg.c);
                }
            }
        }
    }
    // transfer stress from device to host
    syncmem_var_d2h_op()(sigma_onsite.data(), stress_device, 9);
    delmem_var_op()(stress_device);
    // sum up forcenl from all processors
    for (int l = 0; l < 3; l++)
    {
        for (int m = 0; m < 3; m++)
        {
            if (m > l)
            {
                sigma_onsite[l * 3 + m] = sigma_onsite[m * 3 + l];
            }
            Parallel_Reduce::reduce_all(sigma_onsite[l * 3 + m]); // qianrui fix a bug for kpar > 1
        }
    }
    // rescale the stress with 1/omega
    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            sigma_onsite[ipol * 3 + jpol] *= 1.0 / ucell_in.omega;
        }
    }

    for (int ipol = 0; ipol < 3; ipol++)
    {
        for (int jpol = 0; jpol < 3; jpol++)
        {
            sigma(ipol, jpol) = sigma_onsite[ipol * 3 + jpol];
        }
    }
    // do symmetry
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigma, ucell_in.lat);
    } // end symmetry

    ModuleBase::timer::tick("Stress_Func", "stress_onsite");
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif