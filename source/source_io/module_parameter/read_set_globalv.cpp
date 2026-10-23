#include "read_input.h"
#include "read_input_tool.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"

namespace
{
/// @brief Find the structure file for restarting a molecular dynamics run.
/// @param stru_dir directory where the previous run wrote its structure files
/// @param istep number of completed ionic steps recorded in Restart_md.txt
/// @return path of the first existing structure file, empty if none exists
/// @details Two naming schemes are accepted, in order:
/// - "STRU{istep + 1}": written by run_md since the ionic-step numbering
///   was unified with relax (1-based);
/// - "STRU_MD_{istep}": written by older ABACUS versions.
std::string find_md_restart_stru(const std::string& stru_dir, const int istep)
{
    const std::string stru_new = stru_dir + "STRU" + std::to_string(istep + 1);
    std::ifstream file_new(stru_new.c_str());
    if (file_new.good())
    {
        return stru_new;
    }
    const std::string stru_old = stru_dir + "STRU_MD_" + std::to_string(istep);
    std::ifstream file_old(stru_old.c_str());
    if (file_old.good())
    {
        return stru_old;
    }
    return "";
}
} // namespace

namespace ModuleIO
{
/// @note Here para.inp has been synchronized of all ranks.
///       All para.inp have the same value.
void ReadInput::set_globalv(const Input_para& inp, System_para& sys)
{
    /// caculate the gamma_only_pw and gamma_only_local
    if (inp.gamma_only)
    {
        sys.gamma_only_local = true;
    }
    if (sys.gamma_only_local)
    {
        if (inp.esolver_type == "tddft")
        {
            GlobalV::ofs_running << " WARNING : gamma_only is not applicable for tddft" << std::endl;
            sys.gamma_only_local = false;
        }
    }
    /// set deepks_setorb
    if (inp.deepks_scf || inp.deepks_out_labels == 1)
    {
        sys.deepks_setorb = true;
    }
    /// set the noncolin and lspinorb from nspin
    switch (inp.nspin)
    {
    case 4:
        if (inp.noncolin)
        {
            sys.domag = true;
            sys.domag_z = false;
        }
        else
        {
            sys.domag = false;
            sys.domag_z = true;
        }
        sys.npol = 2;
        break;
    case 2:
    case 1:
        sys.domag = false;
        sys.domag_z = false;
        sys.npol = 1;
    default:
        break;
    }
    sys.nqx = static_cast<int>((sqrt(inp.ecutwfc) / sys.dq + 4.0) * inp.cell_factor);
    sys.nqxq = static_cast<int>((sqrt(inp.ecutrho) / sys.dq + 4.0) * inp.cell_factor);
    /// set ncx,ncy,ncz
    sys.ncx = inp.nx;
    sys.ncy = inp.ny;
    sys.ncz = inp.nz;
#ifdef __MPI
    Parallel_Common::bcast_bool(sys.double_grid);
#endif
    /// set ks_run
    if (inp.ks_solver != "bpcg" && inp.bndpar > 1)
    {
        sys.all_ks_run = false;
    }
    // set the has_double_data and has_float_data
#ifdef __FLOAT_FFTW
    bool float_cond = inp.cal_cond && inp.esolver_type == "sdft";
#else
    bool float_cond = false;
#endif
    sys.has_double_data = (inp.precision == "double") || (inp.precision == "mixing") || float_cond;
    sys.has_float_data = (inp.precision == "single") || (inp.precision == "mixing") || float_cond;
}

/// @note Here para.inp has not been synchronized of all ranks.
///       Only para.inp in rank 0 is right.
///       So we need to broadcast the results to all ranks.
void ReadInput::set_global_dir(const Input_para& inp, System_para& sys)
{
    /// caculate the global output directory
    const std::string prefix = "OUT.";
    sys.global_out_dir = prefix + inp.suffix + "/";
    sys.global_out_dir = to_dir(sys.global_out_dir);

    /// Set subdirectories based on out_app_flag
    /// out_app_flag == 2: use unified uppercase subdirectories
    /// otherwise: keep backward-compatible behavior (all point to global_out_dir)
    if (inp.out_app_flag == 2)
    {
        sys.global_stru_dir = sys.global_out_dir + "STRU/";
        sys.global_psi_dir = sys.global_out_dir + "PSI/";
        sys.global_chg_dir = sys.global_out_dir + "CHG/";
        sys.global_pot_dir = sys.global_out_dir + "POT/";
        sys.global_elf_dir = sys.global_out_dir + "ELF/";
        sys.global_dos_dir = sys.global_out_dir + "DOS/";
        sys.global_dip_dir = sys.global_out_dir + "DIP/";
        sys.global_mat_dir = sys.global_out_dir + "MAT/";
    }
    else
    {
        sys.global_stru_dir = sys.global_out_dir;
        sys.global_psi_dir = sys.global_out_dir;
        sys.global_chg_dir = sys.global_out_dir;
        sys.global_pot_dir = sys.global_out_dir;
        sys.global_elf_dir = sys.global_out_dir;
        sys.global_dos_dir = sys.global_out_dir;
        sys.global_dip_dir = sys.global_out_dir;
        sys.global_mat_dir = sys.global_out_dir;
    }

    /// Map legacy directory names to new unified names
    /// WFC -> PSI, matrix -> MAT
    sys.global_wfc_dir = sys.global_psi_dir;
    sys.global_matrix_dir = sys.global_mat_dir;

    /// get the global ML KEDF descriptor directory
    sys.global_mlkedf_descriptor_dir = sys.global_out_dir + "MLKEDF_Descriptors/";
    sys.global_mlkedf_descriptor_dir = to_dir(sys.global_mlkedf_descriptor_dir);

    /// get the global directory for DeePKS labels during electronic steps
    sys.global_deepks_label_elec_dir = sys.global_out_dir + "DeePKS_Labels_Elec/";
    sys.global_deepks_label_elec_dir = to_dir(sys.global_deepks_label_elec_dir);

    /// get the global readin directory
    sys.global_readin_dir = inp.read_file_dir;
    sys.global_readin_dir = to_dir(sys.global_readin_dir);

    /// get the stru file for md restart case
    if (inp.calculation == "md" && inp.mdp.md_restart)
    {
        const int istep = current_md_step(sys.global_readin_dir);

        if (inp.read_file_dir == to_dir("OUT." + inp.suffix))
        {
            // Restart from the default output directory: search the current
            // STRU output directory, which follows out_app_flag.
            const std::string stru_dir = sys.global_stru_dir;
            sys.global_in_stru = find_md_restart_stru(stru_dir, istep);
            if (sys.global_in_stru.empty())
            {
                // Releases before the out_app_flag rework always wrote
                // STRU_MD_{istep} into a fixed STRU/ subdirectory of the
                // output directory; keep those outputs restartable.
                const std::string legacy_dir = sys.global_out_dir + "STRU/";
                sys.global_in_stru = find_md_restart_stru(legacy_dir, istep);
            }
        }
        else
        {
            // Restart from a user-provided directory: the previous run may
            // have written either naming scheme there.
            const std::string stru_dir = inp.read_file_dir;
            sys.global_in_stru = find_md_restart_stru(stru_dir, istep);
        }

        if (sys.global_in_stru.empty())
        {
            const std::string message = "no STRU file found for MD restart in " + sys.global_readin_dir
                                        + "; expected STRU_MD_{step} or STRU{step+1}";
            ModuleBase::WARNING_QUIT("ReadInput::set_global_dir", message);
        }
    }
    else
    {
        sys.global_in_stru = inp.stru_file;
    }

    // set the global log file
    bool out_alllog = inp.out_alllog;
    // set the global calculation type
    std::string cal_type = inp.calculation;
#ifdef __MPI
    // because log_file is different for each rank, so we need to bcast the out_alllog
    Parallel_Common::bcast_bool(out_alllog);
    // In `ReadInput::read_parameters`, `bcastfunc(param)` is after `set_global_dir`,
    // so `cal_type` must be synchronized here manually
    Parallel_Common::bcast_string(cal_type);
#endif
    if (out_alllog)
    {
        PARAM.sys.log_file = "running_" + cal_type + "_" + std::to_string(PARAM.sys.myrank + 1) + ".log";
    }
    else
    {
        PARAM.sys.log_file = "running_" + cal_type + ".log";
    }
#ifdef __MPI
    Parallel_Common::bcast_string(sys.global_in_card);
    Parallel_Common::bcast_string(sys.global_out_dir);
    Parallel_Common::bcast_string(sys.global_readin_dir);
    Parallel_Common::bcast_string(sys.global_stru_dir);
    Parallel_Common::bcast_string(sys.global_matrix_dir);
    Parallel_Common::bcast_string(sys.global_wfc_dir);
    Parallel_Common::bcast_string(sys.global_psi_dir);
    Parallel_Common::bcast_string(sys.global_chg_dir);
    Parallel_Common::bcast_string(sys.global_pot_dir);
    Parallel_Common::bcast_string(sys.global_elf_dir);
    Parallel_Common::bcast_string(sys.global_dos_dir);
    Parallel_Common::bcast_string(sys.global_dip_dir);
    Parallel_Common::bcast_string(sys.global_mat_dir);
    Parallel_Common::bcast_string(sys.global_in_stru);
#endif
}
} // namespace ModuleIO
