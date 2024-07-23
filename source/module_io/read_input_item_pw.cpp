#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_pw()
{
    {
        Input_Item item("ecutwfc");
        item.annotation = "energy cutoff for wave functions";
        read_sync_double(input.ecutwfc);
        this->add_item(item);
    }
    {
        Input_Item item("ecutrho");
        item.annotation = "energy cutoff for charge density and potential";
        read_sync_double(input.ecutrho);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            Input_para& input = para.input;
            if (input.ecutrho <= 0.0)
            {
                input.ecutrho = 4.0 * input.ecutwfc;
            }
            if (input.nx * input.ny * input.nz == 0 && input.ecutrho / input.ecutwfc > 4 + 1e-8)
            {
                para.sys.double_grid = true;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ecutrho / para.input.ecutwfc < 4 - 1e-8)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutrho/ecutwfc must >= 4");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("erf_ecut");
        item.annotation = "the value of the constant energy cutoff";
        read_sync_double(input.erf_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("erf_height");
        item.annotation = "the height of the energy step for reciprocal vectors";
        read_sync_double(input.erf_height);
        this->add_item(item);
    }
    {
        Input_Item item("erf_sigma");
        item.annotation = "the width of the energy step for reciprocal vectors";
        read_sync_double(input.erf_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("fft_mode");
        item.annotation = "mode of FFTW";
        read_sync_int(input.fft_mode);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_nmax");
        item.annotation = "max iteration number for cg";
        read_sync_int(input.pw_diag_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("diago_cg_prec");
        item.annotation = "diago_cg_prec";
        read_sync_int(input.diago_cg_prec);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_ndim");
        item.annotation = "dimension of workspace for Davidson diagonalization";
        read_sync_int(input.pw_diag_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("diago_full_acc");
        item.annotation = "all the empty states are diagonalized";
        read_sync_bool(input.diago_full_acc);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_thr");
        item.annotation = "threshold for eigenvalues is cg electron iterations";
        read_sync_double(input.pw_diag_thr);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_S" && para.input.basis_type == "pw")
            {
                if (para.input.pw_diag_thr > 1.0e-3) {
                    para.input.pw_diag_thr = 1.0e-5;
}
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nb2d");
        item.annotation = "matrix 2d division";
        read_sync_int(input.nb2d);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nb2d < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nb2d should be greater than 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr");
        item.annotation = "charge density error";
        read_sync_double(input.scf_thr);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.scf_thr == -1.0)
            {
                if (para.input.basis_type == "lcao" || para.input.basis_type == "lcao_in_pw")
                {
                    para.input.scf_thr = 1.0e-7;
                }
                else if (para.input.basis_type == "pw" && para.input.calculation != "nscf")
                {
                    para.input.scf_thr = 1.0e-9;
                }
                else if (para.input.basis_type == "pw" && para.input.calculation == "nscf")
                {
                    para.input.scf_thr = 1.0e-6;
                    // In NSCF calculation, the diagonalization threshold is set
                    // to 0.1*scf/nelec. In other words, the scf_thr is used to
                    // control diagonalization convergence threthod in NSCF. In
                    // this case, the default 1.0e-9 is too strict. renxi
                    // 20230908
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr_type");
        item.annotation = "type of the criterion of scf_thr, 1: reci drho for "
                          "pw, 2: real drho for lcao";
        read_sync_int(input.scf_thr_type);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.scf_thr_type == -1)
            {
                if (para.input.basis_type == "lcao" || para.input.basis_type == "lcao_in_pw")
                {
                    para.input.scf_thr_type = 2;
                }
                else if (para.input.basis_type == "pw")
                {
                    para.input.scf_thr_type = 1;
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("init_wfc");
        item.annotation = "start wave functions are from 'atomic', "
                          "'atomic+random', 'random' or";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.init_wfc = "file";
            }
            if (para.input.basis_type == "lcao_in_pw")
            {
                if (para.input.init_wfc != "nao")
                {
                    para.input.init_wfc = "nao";
                    GlobalV::ofs_warning << "init_wfc is set to nao when "
                                            "basis_type is lcao_in_pw"
                                         << std::endl;
                }
            }
        };
        read_sync_string(input.init_wfc);
        this->add_item(item);
    }
    {
        Input_Item item("psi_initializer");
        item.annotation = "whether to use psi_initializer";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.basis_type == "lcao_in_pw")
            {
                para.input.psi_initializer = true;
            }
        };
        read_sync_bool(input.psi_initializer);
        this->add_item(item);
    }
    {
        Input_Item item("init_chg");
        item.annotation = "start charge is from 'atomic' or file";
        read_sync_string(input.init_chg);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.init_chg = "atomic";
            }
            if (para.input.calculation == "nscf" || para.input.calculation == "get_S")
            {
                if (para.input.init_chg != "file")
                {
                    ModuleBase::GlobalFunc::AUTO_SET("init_chg", para.input.init_chg);
                }
                para.input.init_chg = "file";
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.init_chg != "atomic" && para.input.init_chg != "file" && para.input.init_chg != "auto")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "init_chg should be 'atomic', 'file' or 'auto'");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("chg_extrap");
        item.annotation = "atomic; first-order; second-order; dm:coefficients of SIA";
        read_sync_string(input.chg_extrap);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.chg_extrap == "default" && para.input.calculation == "md")
            {
                para.input.chg_extrap = "second-order";
            }
            else if (para.input.chg_extrap == "default"
                     && (para.input.calculation == "relax" || para.input.calculation == "cell-relax"))
            {
                para.input.chg_extrap = "first-order";
            }
            else if (para.input.chg_extrap == "default")
            {
                para.input.chg_extrap = "atomic";
            }
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
                para.input.chg_extrap = "atomic";
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_chg");
        item.annotation = ">0 output charge density for selected electron steps";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg") {
                para.input.out_chg = 1;
}
        };
        read_sync_int(input.out_chg);
        this->add_item(item);
    }
    {
        Input_Item item("out_pot");
        item.annotation = "output realspace potential";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg") {
                para.input.out_pot = 0;
}
        };
        read_sync_int(input.out_pot);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_pw");
        item.annotation = "output wave functions";
        read_sync_int(input.out_wfc_pw);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_r");
        item.annotation = "output wave functions in realspace";
        read_sync_bool(input.out_wfc_r);
        this->add_item(item);
    }
    {
        Input_Item item("out_dos");
        item.annotation = "output energy and dos";
        read_sync_int(input.out_dos);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg") {
                para.input.out_dos = 0;
}
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_dos == 3 && para.input.symmetry == "1")
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "symmetry can't be used for out_dos==3(Fermi Surface "
                                         "Plotting) by now.");
            }
            if (para.input.basis_type == "pw" && para.input.out_dos == 3)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "Fermi Surface Plotting not "
                                         "implemented for plane wave now.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_band");
        item.annotation = "output energy and band structure (with precision 8)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            if (count == 1)
            {
                para.input.out_band[0] = std::stoi(item.str_values[0]);
                para.input.out_band[1] = 8;
            }
            else if (count == 2)
            {
                para.input.out_band[0] = std::stoi(item.str_values[0]);
                para.input.out_band[1] = std::stoi(item.str_values[1]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_band should have 1 or 2 values");
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg") {
                para.input.out_band[0] = 0;
}
        };
        sync_intvec(input.out_band, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_proj_band");
        item.annotation = "output projected band structure";
        read_sync_bool(input.out_proj_band);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg") {
                para.input.out_proj_band = false;
}
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type == "pw" && para.input.out_proj_band)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_proj_band is only for lcao");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("restart_save");
        item.annotation = "print to disk every step for restart";
        read_sync_bool(input.restart_save);
        this->add_item(item);
    }
    {
        Input_Item item("restart_load");
        item.annotation = "restart from disk";
        read_sync_bool(input.restart_load);
        this->add_item(item);
    }
    {
        Input_Item item("read_file_dir");
        item.annotation = "directory of files for reading";
        read_sync_string(input.read_file_dir);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.read_file_dir == "auto")
            {
                para.input.read_file_dir = "OUT." + para.input.suffix;
            }
            else
            {
                para.input.read_file_dir = para.input.read_file_dir;
            }
            para.input.read_file_dir = to_dir(para.input.read_file_dir);
        };
        this->add_item(item);
    }
    {
        Input_Item item("nx");
        item.annotation = "number of points along x axis for FFT grid";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.nx = intvalue;
            para.sys.ncx = intvalue;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.nx != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(input.nx);
        this->add_item(item);
    }
    {
        Input_Item item("ny");
        item.annotation = "number of points along y axis for FFT grid";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.ny = intvalue;
            para.sys.ncy = intvalue;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.ny != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(input.ny);
        this->add_item(item);
    }
    {
        Input_Item item("nz");
        item.annotation = "number of points along z axis for FFT grid";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.nz = intvalue;
            para.sys.ncz = intvalue;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.nz != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(input.nz);
        this->add_item(item);
    }
    {
        Input_Item item("ndx");
        item.annotation = "number of points along x axis for FFT smooth grid";
        read_sync_int(input.ndx);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndx > para.input.nx)
            {
                para.sys.double_grid = true;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            if (para.input.ndx * para.input.ndy * para.input.ndz == 0 && para.input.ndx != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx, ndy, ndz should be all set to non-zero");
            }
            if (para.input.ndx < para.input.nx)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx should be greater than or equal to nx");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("ndy");
        item.annotation = "number of points along y axis for FFT smooth grid";
        read_sync_int(input.ndy);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndy > para.input.ny)
            {
                para.sys.double_grid = true;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            if (para.input.ndx * para.input.ndy * para.input.ndz == 0 && para.input.ndy != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx, ndy, ndz should be all set to non-zero");
            }
            if (para.input.ndy < para.input.ny)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndy should be greater than or equal to ny");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("ndz");
        item.annotation = "number of points along z axis for FFT smooth grid";
        read_sync_int(input.ndz);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndy > para.input.ny)
            {
                para.sys.double_grid = true;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (!item.is_read()) {
                return;
}
            if (para.input.ndx * para.input.ndy * para.input.ndz == 0 && para.input.ndz != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx, ndy, ndz should be all set to non-zero");
            }
            if (para.input.ndz < para.input.nz)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndz should be greater than or equal to nz");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("cell_factor");
        item.annotation = "used in the construction of the pseudopotential tables";
        read_sync_double(input.cell_factor);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "cell-relax" && para.input.cell_factor < 2.0)
            {
                para.input.cell_factor = 2.0; // follows QE
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pw_seed");
        item.annotation = "random seed for initializing wave functions";
        read_sync_int(input.pw_seed);
        this->add_item(item);
    }
}
} // namespace ModuleIO
