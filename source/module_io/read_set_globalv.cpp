#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"
#include "read_input.h"
#include "read_input_tool.h"
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
    if (inp.deepks_scf || inp.deepks_out_labels)
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
}

/// @note Here para.inp has been synchronized of all ranks. 
///       Only para.inp in rank 0 is right. 
///       So we need to broadcast the results to all ranks.
void ReadInput::set_global_dir(const Input_para& inp, System_para& sys)
{
    /// caculate the global output directory
    const std::string prefix = "OUT.";
    sys.global_out_dir = prefix + inp.suffix + "/";
    sys.global_out_dir = to_dir(sys.global_out_dir);

    /// get the global output directory
    sys.global_stru_dir = sys.global_out_dir + "STRU/";
    sys.global_stru_dir = to_dir(sys.global_stru_dir);

    /// get the global output directory
    sys.global_matrix_dir = sys.global_out_dir + "matrix/";
    sys.global_matrix_dir = to_dir(sys.global_matrix_dir);

    /// get the global readin directory
    sys.global_readin_dir = inp.read_file_dir;
    sys.global_readin_dir = to_dir(sys.global_readin_dir);

    /// get the stru file for md restart case
    if (inp.calculation == "md" && inp.mdp.md_restart)
    {
        int istep = current_md_step(sys.global_readin_dir);

        if (inp.read_file_dir == to_dir("OUT." + inp.suffix))
        {
            sys.global_in_stru = sys.global_stru_dir + "STRU_MD_" + std::to_string(istep);
        }
        else
        {
            sys.global_in_stru = inp.read_file_dir + "STRU_MD_" + std::to_string(istep);
        }
    }
    else
    {
        sys.global_in_stru = inp.stru_file;
    }

    // set the global log file
    bool out_alllog = inp.out_alllog;
#ifdef __MPI
    // because log_file is different for each rank, so we need to bcast the out_alllog
    Parallel_Common::bcast_bool(out_alllog);
#endif
    if (out_alllog)
    {
        PARAM.sys.log_file = "running_" + PARAM.inp.calculation + "_" + std::to_string(PARAM.sys.myrank + 1) + ".log";
    }
    else
    {
        PARAM.sys.log_file = "running_" + PARAM.inp.calculation + ".log";
    }
#ifdef __MPI
    Parallel_Common::bcast_string(sys.global_in_card);
    Parallel_Common::bcast_string(sys.global_out_dir);
    Parallel_Common::bcast_string(sys.global_readin_dir);
    Parallel_Common::bcast_string(sys.global_stru_dir);
    Parallel_Common::bcast_string(sys.global_matrix_dir);
    Parallel_Common::bcast_string(sys.global_in_stru);
#endif
}
} // namespace ModuleIO
