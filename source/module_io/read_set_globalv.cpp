#include "module_base/global_variable.h"
#include "module_base/module_device/device.h"
#include "module_base/tool_quit.h"
#include "module_parameter/parameter.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
/// @note Here para.inp has been synchronized of all ranks. 
///       Only para.inp in rank 0 is right. 
///       So we need to broadcast the results to all ranks.
void ReadInput::set_global_dir(Parameter& para)
{
    /// caculate the global output directory
    const std::string prefix = "OUT.";
    para.sys.global_out_dir = prefix + para.inp.suffix + "/";
    para.sys.global_out_dir = to_dir(para.sys.global_out_dir);

    /// get the global output directory
    para.sys.global_stru_dir = para.globalv.global_out_dir + "STRU/";
    para.sys.global_stru_dir = to_dir(para.sys.global_stru_dir);

    /// get the global output directory
    para.sys.global_matrix_dir = para.globalv.global_out_dir + "matrix/";
    para.sys.global_matrix_dir = to_dir(para.sys.global_matrix_dir);

    /// get the global readin directory
    para.sys.global_readin_dir = para.inp.read_file_dir + '/';
    para.sys.global_readin_dir = to_dir(para.sys.global_readin_dir);

    /// get the stru file for md restart case
    if (para.inp.calculation == "md" && para.mdp.md_restart)
    {
        int istep = current_md_step(para.sys.global_readin_dir);

        if (para.inp.read_file_dir == to_dir("OUT." + para.inp.suffix))
        {
            para.sys.global_in_stru = para.sys.global_stru_dir + "STRU_MD_" + std::to_string(istep);
        }
        else
        {
            para.sys.global_in_stru = para.inp.read_file_dir + "STRU_MD_" + std::to_string(istep);
        }
    }
    else
    {
        para.sys.global_in_stru = para.inp.stru_file;
    }

    // set the global log file
    bool out_alllog = para.inp.out_alllog;
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
    Parallel_Common::bcast_string(para.sys.global_in_card);
    Parallel_Common::bcast_string(para.sys.global_out_dir);
    Parallel_Common::bcast_string(para.sys.global_readin_dir);
    Parallel_Common::bcast_string(para.sys.global_stru_dir);
    Parallel_Common::bcast_string(para.sys.global_matrix_dir);
    Parallel_Common::bcast_string(para.sys.global_in_stru);
#endif
}

/// @note Here para.inp has been synchronized of all ranks.
///       All para.inp have the same value.
void ReadInput::set_globalv(Parameter& para)
{
    /// caculate the gamma_only_pw and gamma_only_local
    if (para.inp.gamma_only)
    {
        para.sys.gamma_only_local = true;
    }
    if (para.sys.gamma_only_local)
    {
        if (para.inp.esolver_type == "tddft")
        {
            GlobalV::ofs_running << " WARNING : gamma_only is not applicable for tddft" << std::endl;
            para.sys.gamma_only_local = false;
        }
    }
    /// set deepks_setorb
    if (para.inp.deepks_scf || para.inp.deepks_out_labels)
    {
        para.sys.deepks_setorb = true;
    }
    /// set the noncolin and lspinorb from nspin
    switch (para.inp.nspin)
    {
    case 4:
        if (para.inp.noncolin)
        {
            para.sys.domag = true;
            para.sys.domag_z = false;
        }
        else
        {
            para.sys.domag = false;
            para.sys.domag_z = true;
        }
        para.sys.npol = 2;
        break;
    case 2:
    case 1:
        para.sys.domag = false;
        para.sys.domag_z = false;
        para.sys.npol = 1;
    default:
        break;
    }
    /// set ncx,ncy,ncz
    para.sys.ncx = para.inp.nx;
    para.sys.ncy = para.inp.ny;
    para.sys.ncz = para.inp.nz;
#ifdef __MPI
    Parallel_Common::bcast_bool(para.sys.double_grid);
#endif
    // calculate the number of nbands_local
    para.sys.nbands_l = para.inp.nbands / para.inp.bndpar;
    if (GlobalV::RANK_IN_BPGROUP < para.inp.nbands % para.inp.bndpar)
    {
        para.sys.nbands_l++;
    }
}
} // namespace ModuleIO
