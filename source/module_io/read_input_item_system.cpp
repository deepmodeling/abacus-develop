#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
#include "module_base/module_device/device.h"

#include <fstream>
#include <unistd.h>

namespace ModuleIO
{
// There are some examples:
// Generallly:
// {
//      Input_Item item("suffix");
//      item.annotation = "the name of main output directory";
//      read_sync_string(input.suffix);
//      this->add_item(item);
// }
//
// Specially:
// {
//      Input_Item item("kspacing");
//      item.annotation = "unit in 1/bohr, should be > 0, default is 0 which
//      means read KPT file";
//
//      item.read_value = [](const Input_Item& item, Parameter& para) {
//          para.input.kspacing[0] = std::stod(item.str_values[0]);
//          para.input.kspacing[1] = std::stod(item.str_values[1]);
//          para.input.kspacing[2] = std::stod(item.str_values[2]);
//      };
//
//      item.reset_value = [](const Input_Item& item, Parameter& para) {
//          if(para.input.kspacing[0] <= 0) para.input.kspacing[0] = 1;
//      };
//
//      item.check_value = [](const Input_Item& item, const Parameter& para)
//      {assert(para.input.kspacing[0]>0);};
//
//      item.get_final_value = [](Input_Item& item, const Parameter& para) {
//          item.final_value << para.input.kspacing[0] << " " <<
//          para.input.kspacing[1] << " " << para.input.kspacing[2];
//      };
//
//      add_doublevec_bcast(&Parameter::PARAMETER, N);
//      this->add_item(item);
//  }
void ReadInput::item_system()
{
    {
        Input_Item item("suffix");
        item.annotation = "the name of main output directory";
        read_sync_string(input.suffix);
        this->add_item(item);
    }
    {
        Input_Item item("ntype");
        item.annotation = "atom species number";
        // check of ntype is done in check_ntype
        read_sync_int(input.ntype);
        this->add_item(item);
    }
    {
        Input_Item item("calculation");
        item.annotation = "test; scf; relax; nscf; get_wf; get_pchg";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.calculation = strvalue;
            std::string& calculation = para.input.calculation;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::string& calculation = para.input.calculation;
            std::vector<std::string> callist = {"scf",
                                                "relax",
                                                "md",
                                                "cell-relax",
                                                "test_memory",
                                                "test_neighbour",
                                                "nscf",
                                                "get_S",
                                                "get_wf",
                                                "get_pchg",
                                                "gen_bessel"};
            if (std::find(callist.begin(), callist.end(), calculation) == callist.end())
            {
                const std::string warningstr = nofound_str(callist, "calculation");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
            if (calculation == "get_pchg" || calculation == "get_wf")
            {
                if (para.input.basis_type == "pw") // xiaohui add 2013-09-01
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "calculate = get_pchg or get_wf "
                                             "is only availble for LCAO.");
                }
            }
            else if (calculation == "gen_bessel")
            {
                if (para.input.basis_type != "pw")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "to generate descriptors, please use pw basis");
                }
            }
        };
        sync_string(input.calculation);
        this->add_item(item);
    }
    {
        Input_Item item("esolver_type");
        item.annotation = "the energy solver: ksdft, sdft, ofdft, tddft, lj, dp, ks-lr, lr";
        read_sync_string(input.esolver_type);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> esolver_types = { "ksdft", "sdft", "ofdft", "tddft", "lj", "dp", "lr", "ks-lr" };
            if (std::find(esolver_types.begin(), esolver_types.end(), para.input.esolver_type) == esolver_types.end())
            {
                const std::string warningstr = nofound_str(esolver_types, "esolver_type");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
            if (para.input.esolver_type == "dp")
            {
                if (access(para.input.mdp.pot_file.c_str(), 0) == -1)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "Can not find DP model !");
                }
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.esolver_type == "lr" && para.input.calculation == "scf")
            {   // for LR-only calculation based on the ground-state, set calculation to "nscf"
                para.input.calculation = "nscf";
            }
            };
        this->add_item(item);
    }
    {
        Input_Item item("symmetry");
        item.annotation = "the control of symmetry";
        read_sync_string(input.symmetry);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.symmetry == "default")
            {
                if (para.input.gamma_only || para.input.calculation == "nscf" || para.input.calculation == "get_S"
                    || para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
                {
                    para.input.symmetry = "0"; // if md or exx, symmetry will be
                                               // force-set to 0 or -1 later
                }
                else
                {
                    para.input.symmetry = "1";
                }
            }
            if (para.input.calculation == "md")
            {
                para.input.symmetry = "0";
            }
            if (para.input.efield_flag)
            {
                para.input.symmetry = "0";
            }
            if (para.input.qo_switch)
            {
                para.input.symmetry = "-1"; // disable kpoint reduce
            }
            if (para.input.berry_phase)
            {
                para.input.symmetry = "-1"; // disable kpoint reduce
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_prec");
        item.annotation = "accuracy for symmetry";
        read_sync_double(input.symmetry_prec);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_autoclose");
        item.annotation = "whether to close symmetry automatically when error "
                          "occurs in symmetry analysis";
        read_sync_bool(input.symmetry_autoclose);
        this->add_item(item);
    }
    {
        Input_Item item("cal_stress");
        item.annotation = "calculate the stress or not";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "md")
            {
                if (para.input.esolver_type == "lj" || para.input.esolver_type == "dp"
                    || para.input.mdp.md_type == "msst" || para.input.mdp.md_type == "npt")
                {
                    para.input.cal_stress = true;
                }
            }
            else if (para.input.calculation == "cell-relax")
            {
                para.input.cal_stress = true;
            }
        };
        read_sync_bool(input.cal_stress);
        this->add_item(item);
    }
    {
        Input_Item item("cal_force");
        item.annotation = "if calculate the force at the end of the electronic iteration";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            std::vector<std::string> use_force = {"cell-relax", "relax", "md"};
            std::vector<std::string> not_use_force = {"get_wf", "get_pchg", "nscf", "get_S"};
            if (std::find(use_force.begin(), use_force.end(), para.input.calculation) != use_force.end())
            {
                if (!para.input.cal_force)
                {
                    ModuleBase::GlobalFunc::AUTO_SET("cal_force", "true");
                }
                para.input.cal_force = true;
            }
            else if (std::find(not_use_force.begin(), not_use_force.end(), para.input.calculation) != not_use_force.end())
            {
                if (para.input.cal_force)
                {
                    ModuleBase::GlobalFunc::AUTO_SET("cal_force", "false");
                }
                para.input.cal_force = false;
            }
        };
        read_sync_bool(input.cal_force);
        this->add_item(item);
    }
    {
        Input_Item item("kpar");
        item.annotation = "devide all processors into kpar groups and k points "
                          "will be distributed among";
        read_sync_int(input.kpar);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type == "lcao" && para.input.kpar > 1)
            {
                ModuleBase::WARNING("ReadInput", "kpar > 1 has not been supported for lcao calculation.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bndpar");
        item.annotation = "devide all processors into bndpar groups and bands "
                          "will be distributed among each group";
        read_sync_int(input.bndpar);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.esolver_type != "sdft")
            {
                para.input.bndpar = 1;
            }
            if (para.input.bndpar > GlobalV::NPROC)
            {
                para.input.bndpar = GlobalV::NPROC;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (GlobalV::NPROC % para.input.bndpar != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "The number of processors can not be divided by bndpar");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("latname");
        item.annotation = "the name of lattice name";
        read_sync_string(input.latname);
        this->add_item(item);
    }
    {
        Input_Item item("ecutwfc");
        item.annotation = "energy cutoff for wave functions";
        read_sync_double(input.ecutwfc);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.ecutwfc == 0){ // 0 means no input value
                if (para.input.basis_type == "lcao")
                {
                    para.input.ecutwfc = 100;
                }
                else
                {
                    para.input.ecutwfc = 50;
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ecutwfc <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ecutwfc should be positive");
            }
        };
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
            if (!item.is_read())
            {
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
        Input_Item item("pw_seed");
        item.annotation = "random seed for initializing wave functions";
        read_sync_int(input.pw_seed);
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
            const std::vector<std::string> init_chgs = {"atomic", "file", "wfc", "auto"};
            if (std::find(init_chgs.begin(), init_chgs.end(), para.input.init_chg) == init_chgs.end())
            {
                const std::string warningstr = nofound_str(init_chgs, "init_chg");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dm_to_rho");
        item.annotation = "reads dmr in npz format and calculates electron density";
        read_sync_bool(input.dm_to_rho);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dm_to_rho && GlobalV::NPROC > 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dm_to_rho is not available for parallel calculations");
            }
            if (para.input.dm_to_rho && para.inp.gamma_only)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dm_to_rho is not available for gamma_only calculations");
            }
            if (para.input.dm_to_rho)
            {
#ifndef __USECNPY
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "to write in npz format, please "
                                         "recompile with -DENABLE_CNPY=1");
#endif
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
        Input_Item item("init_vel");
        item.annotation = "read velocity from STRU or not";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "md")
            {
                if (para.input.mdp.md_tfirst < 0 || para.input.mdp.md_restart)
                {
                    para.input.init_vel = true;
                }
            }
        };
        read_sync_bool(input.init_vel);
        this->add_item(item);
    }
    {
        Input_Item item("stru_file");
        item.annotation = "the filename of file containing atom positions";
        read_sync_string(input.stru_file);
        this->add_item(item);
    }
    {
        Input_Item item("kpoint_file");
        item.annotation = "the name of file containing k points";
        read_sync_string(input.kpoint_file);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.stru_file == "")
            {
                GlobalV::ofs_warning << "kpoint_file is set to KPT when stru_file is not set" << std::endl;
                para.input.stru_file = "KPT";
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_dir");
        item.annotation = "the directory containing pseudo files";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            if(item.get_size() == 0)
            {
                para.input.pseudo_dir = "";
            }
            else
            {
                para.input.pseudo_dir = to_dir(strvalue);
            }
        };
        sync_string(input.pseudo_dir);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_dir");
        item.annotation = "the directory containing orbital files";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            if(item.get_size() == 0)
            {
                para.input.orbital_dir = "";
            }
            else
            {
                para.input.orbital_dir = to_dir(strvalue);
            }
        };
        sync_string(input.orbital_dir);
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
            para.input.read_file_dir = to_dir(para.input.read_file_dir);
        };
        this->add_item(item);
    }
    {
        Input_Item item("restart_load");
        item.annotation = "restart from disk";
        read_sync_bool(input.restart_load);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_card");
        item.annotation = "input card for wannier functions";
        read_sync_string(input.wannier_card);
        this->add_item(item);
    }
    {
        Input_Item item("mem_saver");
        item.annotation = "Only for nscf calculations. if set to 1, then a "
                          "memory saving technique will be used for "
                          "many k point calculations.";
        read_sync_int(input.mem_saver);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mem_saver == 1)
            {
                if (para.input.calculation == "scf" || para.input.calculation == "relax")
                {
                    para.input.mem_saver = 0;
                    ModuleBase::GlobalFunc::AUTO_SET("mem_saver", "0");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("diago_proc");
        item.annotation = "the number of procs used to do diagonalization";
        read_sync_int(input.diago_proc);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.diago_proc > GlobalV::NPROC || para.input.diago_proc <= 0)
            {
                para.input.diago_proc = GlobalV::NPROC;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nbspline");
        item.annotation = "the order of B-spline basis";
        read_sync_int(input.nbspline);
        this->add_item(item);
    }
    {
        Input_Item item("kspacing");
        item.annotation = "unit in 1/bohr, should be > 0, default is 0 which "
                          "means read KPT file";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            if (count == 1)
            {
                para.input.kspacing[0] = para.input.kspacing[1] = para.input.kspacing[2] = doublevalue;
            }
            else if (count == 3)
            {
                para.input.kspacing[0] = std::stod(item.str_values[0]);
                para.input.kspacing[1] = std::stod(item.str_values[1]);
                para.input.kspacing[2] = std::stod(item.str_values[2]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "kspacing can only accept one or three values.");
            }
        };
        sync_doublevec(input.kspacing, 3, 0.0);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            int kspacing_zero_num = 0;
            const std::vector<double>& kspacing = para.input.kspacing;
            for (int i = 0; i < 3; i++)
            {
                if (kspacing[i] < 0.0)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "kspacing must > 0");
                }
                else if (kspacing[i] == 0.0)
                {
                    kspacing_zero_num++;
                }
            }
            if (kspacing_zero_num > 0 && kspacing_zero_num < 3)
            {
                std::cout << "kspacing: " << kspacing[0] << " " << kspacing[1] << " " << kspacing[2] << std::endl;
                ModuleBase::WARNING_QUIT("ReadInput", "kspacing must > 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("min_dist_coef");
        item.annotation = "factor related to the allowed minimum distance "
                          "between two atoms";
        read_sync_double(input.min_dist_coef);
        this->add_item(item);
    }
    {
        Input_Item item("device");
        item.annotation = "the computing device for ABACUS";
        read_sync_string(input.device);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            para.input.device=base_device::information::get_device_flag(
                                para.inp.device, para.inp.basis_type);
        };
        this->add_item(item);
    }
    {
        Input_Item item("precision");
        item.annotation = "the computing precision for ABACUS";
        read_sync_string(input.precision);
        this->add_item(item);
    }
}

} // namespace ModuleIO