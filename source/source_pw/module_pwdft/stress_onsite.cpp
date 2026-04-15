#include "source_base/module_device/device.h"
#include "source_base/timer.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "stress_func.h"

/**
 * @brief Calculate the nonlocal pseudopotential stress in PW basis
 * 
 * This function computes the onsite contribution to the stress tensor
 * including DFT+U and spin constraint effects.
 * 
 * @param sigma Stress tensor to be updated
 * @param wg Weight matrix for k-points
 * @param wfc_basis Plane wave basis for wavefunctions
 * @param ucell_in Unit cell information
 * @param dftu DFT+U parameters
 * @param psi_in Wavefunction object (only used for null check)
 * @param p_symm Symmetry object for stress symmetrization
 */
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_onsite(
    ModuleBase::matrix& sigma,
    const ModuleBase::matrix& wg,
    const ModulePW::PW_Basis_K* wfc_basis,
    const UnitCell& ucell_in,
    const Plus_U &dftu, // mohan add 2025-11-06
    const void* psi_in,
    ModuleSymmetry::Symmetry* p_symm
)
{
    ModuleBase::TITLE("Stress", "stress_onsite");
    
    // Early return if required pointers are null
    if (psi_in == nullptr || wfc_basis == nullptr)
    {
        return;
    }
    
    // Check if cell volume is valid (non-zero)
    if (ucell_in.omega <= 0.0)
    {
        ModuleBase::WARNING_QUIT("stress_onsite.cpp", "Cell volume is zero or negative, cannot calculate stress");
    }
    
    ModuleBase::timer::start("Stress", "stress_onsite");

    // Host memory for stress storage (CPU only)
    std::vector<double> sigma_onsite(9, 0.0);

    // Get onsite projector instance
    auto* onsite_projector = projectors::OnsiteProjector<FPTYPE, Device>::get_instance();
    auto* fs_tools = onsite_projector->get_fs_tools();

    const int nks = wfc_basis->nks;
    
    // Loop over all k-points
    for (int ik = 0; ik < nks; ik++)
    {
        // Determine number of occupied bands (skip zero weights)
        int nbands_occ = wg.nc;
        
        // Check if nbands_occ is valid
        if (nbands_occ < 0)
        {
            ModuleBase::WARNING_QUIT("stress_onsite.cpp", "Number of bands is negative, cannot calculate stress");
        }
        
        // Skip if no bands
        if (nbands_occ == 0)
        {
            continue;
        }
        
        // Find the highest occupied band with non-zero weight
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int num_occupied_bands = nbands_occ;

        // Calculate becp = <psi|beta> for all beta functions
        fs_tools->cal_becp(ik, num_occupied_bands);

        // Calculate stress contributions for each tensor component
        for (int ipol = 0; ipol < 3; ipol++)
        {
            for (int jpol = 0; jpol <= ipol; jpol++)
            {
                const int idx = ipol * 3 + jpol;

                // Calculate dbecp_s = <psi|d(beta)/d(epsilon_ij)>
                fs_tools->cal_dbecp_s(ik, num_occupied_bands, ipol, jpol);
                
                // Add DFT+U contribution if enabled
                if (PARAM.inp.dft_plus_u)
                {
                    // Calculate DFT+U stress contribution
                    double dftu_stress = fs_tools->cal_stress_dftu(
                        ik,
                        num_occupied_bands,
                        dftu.orbital_corr.data(),
                        dftu.get_eff_pot_pw(0),
                        dftu.get_size_eff_pot_pw(),
                        wg.c
                    );
                    
                    sigma_onsite[idx] += dftu_stress;
#ifdef __DEBUG
		    std::cout << " idx=" << idx << " stress=" << sigma_onsite[idx] << std::endl;
#endif
                }
                
                // Add spin constraint contribution if enabled
                if (PARAM.inp.sc_mag_switch)
                {
                    // Get spin constraint instance
                    spinconstrain::SpinConstrain<std::complex<double>>& spin_constrain = 
                        spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
                    
                    // Get lambda parameters
                    const std::vector<ModuleBase::Vector3<double>>& lambda = spin_constrain.get_sc_lambda();
                    
                    // Calculate spin constraint stress contribution
                    double dspin_stress = fs_tools->cal_stress_dspin(
                        ik,
                        num_occupied_bands,
                        lambda.data(),
                        wg.c
                    );
                    
                    sigma_onsite[idx] += dspin_stress;
                }
            }
        }
    }

    // No device memory to clean up (using CPU only)

    // Step 1: Reduce stress contributions from all processors
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            const int idx = i * 3 + j;
            if(j>i)
            {
                sigma_onsite[idx]=sigma_onsite[j*3+i];
            }
            Parallel_Reduce::reduce_all(sigma_onsite[idx]); // qianrui fix a bug for kpar > 1
        }
    }

    // Step 2: Rescale stress with 1/omega
    const double inv_omega = 1.0 / ucell_in.omega;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            const int idx = i * 3 + j;
            sigma_onsite[idx] *= inv_omega;
        }
    }

#ifdef __DEBUG
    // Add to total stress
    for (int idx = 0; idx < 9; idx++)
    {
	    std::cout << " idx=" << idx << " stress=" << sigma_onsite[idx] << std::endl;
    }
#endif

    // Step 3: Assign stress values to output matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            const int idx = i * 3 + j;
            sigma(i, j) = sigma_onsite[idx];
        }
    }

    // Symmetrize stress tensor if symmetry is enabled
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        p_symm->symmetrize_mat3(sigma, ucell_in.lat);
    }

    ModuleBase::timer::end("Stress", "stress_onsite");
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif
