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
        read_sync_double(ecutwfc);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.ecutrho <= 0.0)
            {
                para.input.ecutrho = 4.0 * para.input.ecutwfc;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("ecutrho");
        item.annotation = "energy cutoff for charge density and potential";
        read_sync_double(ecutrho);
        autosetfuncs.push_back([](Parameter& para) {
            Input_para& input = para.input;
            if (input.ecutrho <= 0.0)
            {
                input.ecutrho = 4.0 * input.ecutwfc;
            }
            if (input.nx * input.ny * input.nz == 0 && input.ecutrho / input.ecutwfc > 4 + 1e-8)
            {
                input.sup.double_grid = true;
            }
        });
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ecutrho / para.input.ecutwfc < 4 - 1e-8)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutrho/ecutwfc must >= 4");
            }
        };
        add_bool_bcast(sup.double_grid);
        this->add_item(item);
    }
    {
        Input_Item item("erf_ecut");
        item.annotation = "the value of the constant energy cutoff";
        read_sync_double(erf_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("erf_height");
        item.annotation = "the height of the energy step for reciprocal vectors";
        read_sync_double(erf_height);
        this->add_item(item);
    }
    {
        Input_Item item("erf_sigma");
        item.annotation = "the width of the energy step for reciprocal vectors";
        read_sync_double(erf_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("fft_mode");
        item.annotation = "mode of FFTW";
        read_sync_int(fft_mode);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_nmax");
        item.annotation = "max iteration number for cg";
        read_sync_int(pw_diag_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("diago_cg_prec");
        item.annotation = "diago_cg_prec";
        read_sync_int(diago_cg_prec);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_ndim");
        item.annotation = "dimension of workspace for Davidson diagonalization";
        read_sync_int(pw_diag_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("diago_full_acc");
        item.annotation = "all the empty states are diagonalized";
        read_sync_bool(diago_full_acc);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_thr");
        item.annotation = "threshold for eigenvalues is cg electron iterations";
        read_sync_double(pw_diag_thr);
        this->add_item(item);
    }
    {
        Input_Item item("nb2d");
        item.annotation = "matrix 2d division";
        read_sync_int(nb2d);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
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
        read_sync_double(scf_thr);
        autosetfuncs.push_back([](Parameter& para) {
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
                    // In NSCF calculation, the diagonalization threshold is set to 0.1*scf/nelec.
                    // In other words, the scf_thr is used to control diagonalization convergence
                    // threthod in NSCF. In this case, the default 1.0e-9 is too strict.
                    // renxi 20230908
                }
            }
        });
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr_type");
        item.annotation = "type of the criterion of scf_thr, 1: reci drho for pw, 2: real drho for lcao";
        read_sync_int(scf_thr_type);
        autosetfuncs.push_back([](Parameter& para) {
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
        });
        this->add_item(item);
    }
    {
        Input_Item item("init_wfc");
        item.annotation = "start wave functions are from 'atomic', 'atomic+random', 'random' or";
        read_sync_string(init_wfc);
        this->add_item(item);
    }
    {
        Input_Item item("psi_initializer");
        item.annotation = "whether to use psi_initializer";
        read_sync_bool(psi_initializer);
        this->add_item(item);
    }
    {
        Input_Item item("init_chg");
        item.annotation = "start charge is from 'atomic' or file";
        read_sync_string(init_chg);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.init_chg != "atomic" && para.input.init_chg != "file")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "init_chg should be 'atomic' or 'file'");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("chg_extrap");
        item.annotation = "atomic; first-order; second-order; dm:coefficients of SIA";
        read_sync_string(chg_extrap);
        autosetfuncs.push_back([this](Parameter& para) {
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
        });
        this->add_item(item);
    }
    {
        Input_Item item("out_chg");
        item.annotation = ">0 output charge density for selected electron steps";
        read_sync_int(out_chg);
        this->add_item(item);
    }
    {
        Input_Item item("out_pot");
        item.annotation = "output realspace potential";
        read_sync_int(out_pot);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_pw");
        item.annotation = "output wave functions";
        read_sync_int(out_wfc_pw);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_r");
        item.annotation = "output wave functions in realspace";
        read_sync_bool(out_wfc_r);
        this->add_item(item);
    }
    {
        Input_Item item("out_dos");
        item.annotation = "output energy and dos";
        read_sync_int(out_dos);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type == "pw" && para.input.out_dos == 3)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "Fermi Surface Plotting not implemented for plane wave now.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_band");
        item.annotation = "output energy and band structure (with precision 8)";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
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
        sync_intvec(out_band, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_proj_band");
        item.annotation = "output projected band structure";
        read_sync_bool(out_proj_band);
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
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
        read_sync_bool(restart_save);
        this->add_item(item);
    }
    {
        Input_Item item("restart_load");
        item.annotation = "restart from disk";
        read_sync_bool(restart_load);
        this->add_item(item);
    }
    {
        Input_Item item("read_file_dir");
        item.annotation = "directory of files for reading";
        read_sync_string(read_file_dir);
        autosetfuncs.push_back([](Parameter& para) {
            if (para.input.read_file_dir == "auto")
            {
                para.input.sup.readin_dir = "OUT." + para.input.suffix + "/";
            }
            else
            {
                para.input.sup.readin_dir = para.input.read_file_dir + '/';
            }
        });
        this->add_item(item);
    }
    {
        Input_Item item("nx");
        item.annotation = "number of points along x axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.nx = intvalue;
            para.input.sup.ncx = intvalue;
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.nx != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(nx);
        add_int_bcast(sup.ncx);
        this->add_item(item);
    }
    {
        Input_Item item("ny");
        item.annotation = "number of points along y axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.ny = intvalue;
            para.input.sup.ncy = intvalue;
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.ny != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(ny);
        add_int_bcast(sup.ncy);
        this->add_item(item);
    }
    {
        Input_Item item("nz");
        item.annotation = "number of points along z axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.nz = intvalue;
            para.input.sup.ncz = intvalue;
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nx * para.input.ny * para.input.nz == 0 && para.input.nz != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nx, ny, nz should be all set to non-zero");
            }
        };
        sync_int(nz);
        add_int_bcast(sup.ncz);
        this->add_item(item);
    }
    {
        Input_Item item("ndx");
        item.annotation = "number of points along x axis for FFT smooth grid";
        read_sync_int(ndx);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndx > para.input.nx)
            {
                para.input.sup.double_grid = true;
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ndx * para.input.ndy * para.input.ndz == 0 && para.input.ndx != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx, ndy, ndz should be all set to non-zero");
            }
            if (para.input.ndx < para.input.nx)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ndx should be greater than or equal to nx");
            }
        };
        add_bool_bcast(sup.double_grid);
        this->add_item(item);
    }
    {
        Input_Item item("ndy");
        item.annotation = "number of points along y axis for FFT smooth grid";
        read_sync_int(ndy);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndy > para.input.ny)
            {
                para.input.sup.double_grid = true;
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
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
        read_sync_int(ndz);
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if (para.input.ndy > para.input.ny)
            {
                para.input.sup.double_grid = true;
            }
        };
        item.checkvalue = [](const Input_Item& item, const Parameter& para) {
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
        read_sync_double(cell_factor);
        this->add_item(item);
    }
    {
        Input_Item item("pw_seed");
        item.annotation = "random seed for initializing wave functions";
        read_sync_int(pw_seed);
        this->add_item(item);
    }
}
} // namespace ModuleIO
