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
        /// set deepks_setorb
        if (para.input.deepks_scf || para.input.deepks_out_labels)
        {
            para.sys.deepks_setorb = true;
        }
        switch (para.input.nspin)
        {
        case 4:
            if (para.input.noncolin)
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
            para.input.lspinorb = false;
            para.input.noncolin = false;
        default:
            break;
        }
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

    add_bool_bcast(sys.deepks_setorb);
    
    add_bool_bcast(sys.domag);
    add_bool_bcast(sys.domag_z);
    add_int_bcast(sys.npol);

    add_bool_bcast(sys.double_grid);
    add_double_bcast(sys.uramping);
}
} // namespace ModuleIO
