#include "read_input.h"
#include "read_input_tool.h"
#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_base/tool_quit.h"
namespace ModuleIO
{
void ReadInput::set_globalv(Parameter& para)
{
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
        if (para.inp.read_file_dir == "auto")
        {
            para.sys.global_readin_dir = para.globalv.global_out_dir;
        }
        else
        {
            para.sys.global_readin_dir = para.inp.read_file_dir + '/';
        }
        para.sys.global_readin_dir = to_dir(para.sys.global_readin_dir);
        /// caculate the gamma_only_pw and gamma_only_local
        if (para.input.gamma_only)
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

        if ((para.inp.out_mat_r || para.inp.out_mat_hs2 || para.inp.out_mat_t 
                || para.inp.out_mat_dh || para.inp.out_hr_npz
                || para.inp.out_dm_npz || para.inp.dm_to_rho)
            && para.sys.gamma_only_local)
        {
            ModuleBase::WARNING_QUIT("ReadInput",
                                        "output of r(R)/H(R)/S(R)/T(R)/dH(R)/DM(R) is not "
                                        "available for gamma only calculations");
        }

        para.sys.nqx=static_cast<int>((sqrt(para.inp.ecutwfc) / para.sys.dq + 4.0) * para.inp.cell_factor); 
        para.sys.nqxq=static_cast<int>((sqrt(para.inp.ecutrho) / para.sys.dq + 4.0) * para.inp.cell_factor);
    }
}

void ReadInput::set_globalv_bcast()
{
    add_bool_bcast(sys.two_fermi);
    add_bool_bcast(sys.dos_setemin);
    add_bool_bcast(sys.dos_setemax);

    add_int_bcast(sys.ncx);
    add_int_bcast(sys.ncy);
    add_int_bcast(sys.ncz);
    add_bool_bcast(sys.out_md_control);
    add_bool_bcast(sys.rpa_setorb);

    add_bool_bcast(sys.gamma_only_pw);
    add_bool_bcast(sys.gamma_only_local);
    
    add_string_bcast(sys.global_in_card);
    add_string_bcast(sys.global_out_dir);
    add_string_bcast(sys.global_readin_dir);
    add_string_bcast(sys.global_stru_dir);
    add_string_bcast(sys.global_matrix_dir);

    add_bool_bcast(sys.double_grid);
    add_double_bcast(sys.uramping);
    add_double_bcast(sys.dq);
    add_int_bcast(sys.nqx);
    add_int_bcast(sys.nqxq);
}
} // namespace ModuleIO
