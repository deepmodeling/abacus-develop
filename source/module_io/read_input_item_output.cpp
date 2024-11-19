#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_output()
{
    {
        Input_Item item("out_stru");
        item.annotation = "output the structure files after each ion step";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::vector<std::string> offlist = {"nscf", "get_S", "get_pchg", "get_wf"};
            if (std::find(offlist.begin(), offlist.end(), para.input.calculation) != offlist.end())
            {
                para.input.out_stru = false;
            }
        };
        read_sync_bool(input.out_stru);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_elec");
        item.annotation = "the frequency ( >= 0) of electronic iter to output "
                          "charge density and wavefunction. 0: "
                          "output only when converged";
        read_sync_int(input.out_freq_elec);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_ion");
        item.annotation = "the frequency ( >= 0 ) of ionic step to output "
                          "charge density and wavefunction. 0: output "
                          "only when ion steps are finished";
        read_sync_int(input.out_freq_ion);
        this->add_item(item);
    }
    {
        Input_Item item("out_chg");
        item.annotation = "> 0 output charge density for selected electron steps"
                          ", second parameter controls the precision, default is 3.";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_chg should have 1 or 2 values");
            }
            para.input.out_chg[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_chg[1] = (count == 2) ? std::stoi(item.str_values[1]) : 3;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            para.input.out_chg[0] = (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
                                     ? 1
                                     : para.input.out_chg[0];
        };
        sync_intvec(input.out_chg, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_pot");
        item.annotation = "output realspace potential";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
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
        Input_Item item("printe");
        item.annotation = "Print out energy for each band for every printe steps";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.printe <= 0) // default is scf_nmax
            {
                para.input.printe = para.input.scf_nmax;
            }
        };
        read_sync_int(input.printe);
        this->add_item(item);
    }
    {
        Input_Item item("out_band");
        item.annotation = "output energy and band structure (with precision 8)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_band should have 1 or 2 values");
            }
            para.input.out_band[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_band[1] = (count == 2) ? std::stoi(item.str_values[1]) : 8;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
                para.input.out_band[0] = 0;
            }
        };
        sync_intvec(input.out_band, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_dos");
        item.annotation = "output energy and dos";
        read_sync_int(input.out_dos);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
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
        Input_Item item("out_mul");
        item.annotation = "mulliken charge or not";
        read_sync_bool(input.out_mul);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.basis_type == "pw" && para.input.out_mul)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mul is only for lcao");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_proj_band");
        item.annotation = "output projected band structure";
        read_sync_bool(input.out_proj_band);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
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
        Input_Item item("out_level");
        item.annotation = "ie(for electrons); i(for ions);";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.out_level = strvalue;
            para.sys.out_md_control = true;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (!para.sys.out_md_control && para.input.calculation == "md")
            {
                para.input.out_level = "m"; // zhengdy add 2019-04-07
            }
        };
        sync_string(input.out_level);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm");
        item.annotation = ">0 output density matrix";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.out_dm = false;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.sys.gamma_only_local == false && para.input.out_dm)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_dm with k-point algorithm is not implemented yet.");
            }
        };
        read_sync_bool(input.out_dm);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm1");
        item.annotation = ">0 output density matrix (multi-k points)";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_pchg" || para.input.calculation == "get_wf")
            {
                para.input.out_dm1 = false;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.sys.gamma_only_local == true && para.input.out_dm1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_dm1 is only for multi-k");
            }
        };
        read_sync_bool(input.out_dm1);
        this->add_item(item);
    }
    {
        Input_Item item("out_bandgap");
        item.annotation = "if true, print out bandgap";
        read_sync_bool(input.out_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs");
        item.annotation = "output H and S matrix (with precision 8)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_hs should have 1 or 2 values");
            }
            para.input.out_mat_hs[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_mat_hs[1] = (count == 2) ? std::stoi(item.str_values[1]) : 8;
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.qo_switch)
            {
                para.input.out_mat_hs[0] = 1; // print H(k) and S(k)
            }
        };
        sync_intvec(input.out_mat_hs, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_tk");
        item.annotation = "output T(k)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_tk should have 1 or 2 values");
            }
            para.input.out_mat_tk[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_mat_tk[1] = (count == 2) ? std::stoi(item.str_values[1]) : 8;
        };
        sync_intvec(input.out_mat_tk, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs2");
        item.annotation = "output H(R) and S(R) matrix";
        read_sync_bool(input.out_mat_hs2);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_mat_r && para.sys.gamma_only_local)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_r is not available for gamma only calculations");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_dh");
        item.annotation = "output of derivative of H(R) matrix";
        read_sync_bool(input.out_mat_dh);
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
        read_sync_bool(input.out_mat_xc);
        this->add_item(item);
    }
    {
        Input_Item item("out_eband_terms");
        item.annotation = "output the band energy terms separately";
        read_sync_bool(input.out_eband_terms);
        this->add_item(item);
    }
    {
        Input_Item item("out_hr_npz");
        item.annotation = "output hr(I0,JR) submatrices in npz format";
        read_sync_bool(input.out_hr_npz);
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
        read_sync_bool(input.out_dm_npz);
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
        Input_Item item("out_interval");
        item.annotation = "interval for printing H(R) and S(R) matrix during MD";
        read_sync_int(input.out_interval);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_interval <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_interval should be larger than 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_app_flag");
        item.annotation = "whether output r(R), H(R), S(R), T(R), and dH(R) "
                          "matrices in an append manner during MD";
        read_sync_bool(input.out_app_flag);
        this->add_item(item);
    }
    {
        Input_Item item("out_ndigits");
        item.annotation = "the length of decimal part of output data";
        read_sync_int(input.out_ndigits);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_t");
        item.annotation = "output T(R) matrix";
        read_sync_bool(input.out_mat_t);
        this->add_item(item);
    }
    {
        Input_Item item("out_element_info");
        item.annotation = "output (projected) wavefunction of each element";
        read_sync_bool(input.out_element_info);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_r");
        item.annotation = "output r(R) matrix";
        read_sync_bool(input.out_mat_r);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.inp.out_mat_r || para.inp.out_mat_hs2 || para.inp.out_mat_t 
                    || para.inp.out_mat_dh || para.inp.out_hr_npz
                    || para.inp.out_dm_npz || para.inp.dm_to_rho)
                && para.sys.gamma_only_local)
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                            "output of r(R)/H(R)/S(R)/T(R)/dH(R)/DM(R) is not "
                                            "available for gamma only calculations");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lcao");
        item.annotation = "ouput LCAO wave functions, 0, no output 1: text, 2: binary";
        read_sync_int(input.out_wfc_lcao);
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
        Input_Item item("out_dipole");
        item.annotation = "output dipole or not";
        read_sync_bool(input.out_dipole);
        this->add_item(item);
    }
    {
        Input_Item item("out_efield");
        item.annotation = "output dipole or not";
        read_sync_bool(input.out_efield);
        this->add_item(item);
    }
    {
        Input_Item item("out_current");
        item.annotation = "output current or not";
        read_sync_bool(input.out_current);
        this->add_item(item);
    }
    {
        Input_Item item("out_current_k");
        item.annotation = "output current for each k";
        read_sync_bool(input.out_current_k);
        this->add_item(item);
    }
    {
        Input_Item item("out_vecpot");
        item.annotation = "output TDDFT vector potential or not";
        read_sync_bool(input.out_vecpot);
        this->add_item(item);
    }
    {
        Input_Item item("restart_save");
        item.annotation = "print to disk every step for restart";
        read_sync_bool(input.restart_save);
        this->add_item(item);
    }
    {
        Input_Item item("rpa");
        item.annotation = "true:generate output files used in rpa calculation; "
                          "false:(default)";
        read_sync_bool(input.rpa);
        this->add_item(item);
    }
    {
        Input_Item item("nbands_istate");
        item.annotation = "number of bands around Fermi level for get_wf and get_pchg calulation";
        read_sync_int(input.nbands_istate);
        this->add_item(item);
    }
    {
        Input_Item item("bands_to_print");
        item.annotation = "specify the bands to be calculated for the partial (band-decomposed) charge densities";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.bands_to_print);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_intvec_bcast(input.bands_to_print, para.input.bands_to_print.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_pchg");
        item.annotation = "specify the bands to be calculated for the partial (band-decomposed) charge densities";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.out_pchg);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_intvec_bcast(input.out_pchg, para.input.out_pchg.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_norm");
        item.annotation = "specify the bands to be calculated for the norm of wavefunctions";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.out_wfc_norm);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_intvec_bcast(input.out_wfc_norm, para.input.out_wfc_norm.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_re_im");
        item.annotation = "specify the bands to be calculated for the real and imaginary parts of wavefunctions";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.out_wfc_re_im);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_intvec_bcast(input.out_wfc_re_im, para.input.out_wfc_re_im.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("if_separate_k");
        item.annotation = "specify whether to write the partial charge densities for all k-points to individual files "
                          "or merge them";
        read_sync_bool(input.if_separate_k);
        this->add_item(item);
    }
    {
        Input_Item item("out_elf");
        item.annotation = "> 0 output electron localization function (ELF) for selected electron steps"
                          ", second parameter controls the precision, default is 3.";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            std::vector<int> out_elf(count); // create a placeholder vector
            std::transform(item.str_values.begin(), item.str_values.end(), out_elf.begin(), [](std::string s) { return std::stoi(s); });
            // assign non-negative values to para.input.out_elf
            std::copy(out_elf.begin(), out_elf.end(), para.input.out_elf.begin());
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_elf[0] > 0 && para.input.esolver_type != "ksdft" && para.input.esolver_type != "ofdft")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ELF is only aviailable for ksdft and ofdft");
            }
        };
        sync_intvec(input.out_elf, 2, 0);
        this->add_item(item);
    }
}
} // namespace ModuleIO
