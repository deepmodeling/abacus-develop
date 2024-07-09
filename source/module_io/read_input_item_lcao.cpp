
#include <fstream>

#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_lcao()
{
    {
        Input_Item item("basis_type");
        item.annotation = "PW; LCAO in pw; LCAO";
        read_sync_string(basis_type);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.towannier90)
            {
                if (para.input.basis_type == "lcao_in_pw")
                {
                    para.input.basis_type = "lcao";
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> basis_types = {"pw", "lcao_in_pw", "lcao"};
            if (!find_str(basis_types, para.input.basis_type))
            {
                ModuleBase::WARNING_QUIT("ReadInput", "basis_type should be pw, lcao_in_pw, or lcao");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("gamma_only");
        item.annotation = "Only for localized orbitals set and gamma point. If "
                          "set to 1, a fast algorithm is used";
        read_sync_bool(gamma_only);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            Input_para& input = para.input;
            std::string& basis_type = input.basis_type;
            bool& gamma_only = input.gamma_only;
            if (basis_type == "pw" && gamma_only) // pengfei Li add 2015-1-31
            {
                gamma_only = false;
                GlobalV::ofs_warning << " WARNING : gamma_only has not been "
                                        "implemented for pw yet"
                                     << std::endl;
                GlobalV::ofs_warning << " the INPUT parameter gamma_only has been reset to 0" << std::endl;
                GlobalV::ofs_warning << " and a new KPT is generated with "
                                        "gamma point as the only k point"
                                     << std::endl;

                GlobalV::ofs_warning << " Auto generating k-points file: " << input.kpoint_file << std::endl;
                std::ofstream ofs(input.kpoint_file.c_str());
                ofs << "K_POINTS" << std::endl;
                ofs << "0" << std::endl;
                ofs << "Gamma" << std::endl;
                ofs << "1 1 1 0 0 0" << std::endl;
                ofs.close();
            }
            else if (basis_type == "lcao" && gamma_only == 1)
            {
                input.sup.gamma_only_local = true;
                // std::cout << "gamma_only_local =" << gamma_only_local <<
                // std::endl;
                if (input.esolver_type == "tddft")
                {
                    GlobalV::ofs_running << " WARNING : gamma_only is not applicable for tddft" << std::endl;
                    input.sup.gamma_only_local = false;
                }
            }

            if ((input.out_mat_r || input.out_mat_hs2 || input.out_mat_t || input.out_mat_dh || input.out_hr_npz
                 || input.out_dm_npz || input.dm_to_rho)
                && input.sup.gamma_only_local)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "output of r(R)/H(R)/S(R)/T(R)/dH(R)/DM(R) is not "
                                         "available for gamma only calculations");
            }
        };
        add_bool_bcast(sup.gamma_only_local);
        this->add_item(item);
    }
    {
        Input_Item item("search_radius");
        item.annotation = "input search radius (Bohr)";
        read_sync_double(search_radius);
        this->add_item(item);
    }
    {
        Input_Item item("search_pbc");
        item.annotation = "input periodic boundary condition";
        read_sync_bool(search_pbc);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_ecut");
        item.annotation = "energy cutoff for LCAO";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.lcao_ecut == 0 && para.input.basis_type == "lcao")
            {
                para.input.lcao_ecut = para.input.ecutwfc;
                ModuleBase::GlobalFunc::AUTO_SET("lcao_ecut", para.input.ecutwfc);
            }
        };
        read_sync_double(lcao_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dk");
        item.annotation = "delta k for 1D integration in LCAO";
        read_sync_double(lcao_dk);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dr");
        item.annotation = "delta r for 1D integration in LCAO";
        read_sync_double(lcao_dr);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_rmax");
        item.annotation = "max R for 1D two-center integration table";
        read_sync_double(lcao_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs");
        item.annotation = "output H and S matrix (with precision 8)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            if (count == 1)
            {
                para.input.out_mat_hs[0] = std::stoi(item.str_values[0]);
                para.input.out_mat_hs[1] = 8;
            }
            else if (count == 2)
            {
                para.input.out_mat_hs[0] = std::stoi(item.str_values[0]);
                para.input.out_mat_hs[1] = std::stoi(item.str_values[1]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_hs should have 1 or 2 values");
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_switch)
            {
                para.input.out_mat_hs[0] = 1; // print H(k) and S(k)
            }
        };
        sync_intvec(out_mat_hs, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs2");
        item.annotation = "output H(R) and S(R) matrix";
        read_sync_bool(out_mat_hs2);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_dh");
        item.annotation = "output of derivative of H(R) matrix";
        read_sync_bool(out_mat_dh);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_mat_dh && para.input.nspin == 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_dh is not available for nspin = 4");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_xc");
        item.annotation = "output exchange-correlation matrix in KS-orbital representation";
        read_sync_bool(out_mat_xc);
        this->add_item(item);
    }
    {
        Input_Item item("out_hr_npz");
        item.annotation = "output hr(I0,JR) submatrices in npz format";
        read_sync_bool(out_hr_npz);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_hr_npz)
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
        Input_Item item("out_dm_npz");
        item.annotation = "output dmr(I0,JR) submatrices in npz format";
        read_sync_bool(out_dm_npz);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_dm_npz)
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
        Input_Item item("dm_to_rho");
        item.annotation = "reads dmr in npz format and calculates electron density";
        read_sync_bool(dm_to_rho);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dm_to_rho && GlobalV::NPROC > 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dm_to_rho is not available for parallel calculations");
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
        Input_Item item("out_interval");
        item.annotation = "interval for printing H(R) and S(R) matrix during MD";
        read_sync_int(out_interval);
        this->add_item(item);
    }
    {
        Input_Item item("out_app_flag");
        item.annotation = "whether output r(R), H(R), S(R), T(R), and dH(R) "
                          "matrices in an append manner during MD";
        read_sync_bool(out_app_flag);
        this->add_item(item);
    }
    {
        Input_Item item("out_ndigits");
        item.annotation = "the length of decimal part of output data";
        read_sync_int(out_ndigits);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_t");
        item.annotation = "output T(R) matrix";
        read_sync_bool(out_mat_t);
        this->add_item(item);
    }
    {
        Input_Item item("out_element_info");
        item.annotation = "output (projected) wavefunction of each element";
        read_sync_bool(out_element_info);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_r");
        item.annotation = "output r(R) matrix";
        read_sync_bool(out_mat_r);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lcao");
        item.annotation = "ouput LCAO wave functions, 0, no output 1: text, 2: binary";
        read_sync_int(out_wfc_lcao);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_switch)
            {
                para.input.out_wfc_lcao = 1;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_wfc_lcao < 0 || para.input.out_wfc_lcao > 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_wfc_lcao should be 0, 1, or 2");
            }
            if (para.input.basis_type != "lcao" && para.input.out_wfc_lcao != 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_wfc_lcao is only available for basis_type = lcao");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bx");
        item.annotation = "division of an element grid in FFT grid along x";
        read_sync_int(bx);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bx > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bx should be no more than 10");
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.basis_type == "pw" || para.input.basis_type == "lcao_in_pw"
                || para.input.calculation == "get_wf")
            {
                para.input.bx = 1;
                para.input.by = 1;
                para.input.bz = 1;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("by");
        item.annotation = "division of an element grid in FFT grid along y";
        read_sync_int(by);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.by > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "by should be no more than 10");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bz");
        item.annotation = "division of an element grid in FFT grid along z";
        read_sync_int(bz);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bz > 10)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bz should be no more than 10");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("num_stream");
        item.annotation = "the nstream in compute the LCAO with CUDA";
        read_sync_int(nstream);
        this->add_item(item);
    }
    {
        Input_Item item("elpa_num_thread");
        item.annotation = "Number of threads need to use in elpa";
        read_sync_int(elpa_num_thread);
        this->add_item(item);
    }
}
} // namespace ModuleIO