#include "source_pw/module_pwdft/setup_dftu_pw.h"
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 2025-11-06
#include "source_pw/module_pwdft/dftu_base_io.h" // mohan add 2025-11-08
#include "source_pw/module_pwdft/dftu_pw.h"
#include "source_pw/module_pwdft/dftu_pw_tools.h"

namespace DFTU_BASE
{

void iter_init_dftu_pw(const int iter,
                        const int istep,
                        Plus_U_Base& dftu, // mohan add 2025-11-06
                        const void* psi,
                        const ModuleBase::matrix& wg,
                        const UnitCell& ucell,
                        Charge_Mixing* p_chgmix,
                        const std::string& global_out_dir,
                        const OccmatOutputCfg& occmat_cfg,
                        const int* isk,
                        const int nspin,
                        const int kpar)
{
    if (!p_chgmix || occmat_cfg.dft_plus_u <= 0)
    {
        return;
    }

    // Prepare the per-ionic-step file before the first snapshot is appended.
    DFTU_BASE::prepare_ion_step_output(global_out_dir, istep, iter, occmat_cfg);

    const bool first_scf_step = (iter == 1 && istep == 0);

    // Build the +U potential from the file-loaded occupation matrix when:
    //   - init_occ_mat == 2: the matrix is fixed for the whole run, so
    //     rebuild pot_onsite and E_plusU from it at every electronic step
    //     (the matrix is never reprojected from psi or mixed);
    //   - init_occ_mat == 1 at the first SCF step: the file matrix seeds
    //     the Hamiltonian so the first diagonalization includes the +U
    //     correction.  On later steps cal_occ_pw reprojects the matrix
    //     from psi (which now carries the +U influence) and rebuilds the
    //     potential, matching the documented "first step from file, then
    //     update" semantics.  Without this the file matrix would never
    //     enter the Hamiltonian (default mixing_dftu == 0).
    if (dftu.get_init_occ_mat() == 2
        || (first_scf_step && dftu.get_init_occ_mat() == 1))
    {
        pw::compute_pot_uterm_and_energy(ucell,
                                         nspin,
                                         dftu.get_u_current_vec(),
                                         dftu.get_l_channel_vec(),
                                         dftu.get_uterm_mat_index(),
                                         dftu.occmat(),
                                         dftu.get_uterm_mat(),
                                         dftu.energy_ref());
    }
    else if (!first_scf_step)
    {
        // No occupation matrix exists before the first projection, so the
        // matrix can only be accumulated from psi on later electronic steps.
        DFTU_BASE::cal_occ_pw(psi, wg, ucell, p_chgmix, isk, kpar,
                              nspin, dftu.get_device(),
                              dftu.get_l_channel_vec(), dftu.get_u_current_vec(),
                              dftu.get_uterm_mat_index(),
                              dftu.occmat(),
                              dftu.has_occ_mixer() ? &dftu.occ_mixer() : nullptr,
                              dftu.get_uterm_mat(), dftu.energy_ref());
    }

}

}
