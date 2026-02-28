#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_output()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("out_freq_ion");
        item.annotation = "print information every few ionic steps";
        item.category = "Output information";
        item.type = "Integer";
        item.description = "Controls the output interval in ionic steps. When set to a positive integer, information such as charge density, local potential, electrostatic potential, Hamiltonian matrix, overlap matrix, density matrix, and Mulliken population analysis is printed every n ionic steps."
                          "\n\n[NOTE] In RT-TDDFT calculations, this parameter is inactive; output frequency is instead controlled by out_freq_td.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.out_freq_ion <= 0)
            {
                para.input.out_freq_ion = 0; // 0 means no output of info
            }
        };
        read_sync_int(input.out_freq_ion);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_td");
        item.annotation = "print information every few completed electronic iterations in RT-TDDFT";
        item.category = "Output information";
        item.type = "Integer";
        item.description = "Controls the output interval in completed electronic evolution steps during RT-TDDFT calculations. When set to a positive integer n, detailed information (see out_freq_ion) is printed every n electron time-evolution steps (i.e., every STEP OF ELECTRON EVOLVE). For example, if you wish to output information once per ionic step, you should set out_freq_td equal to estep_per_md, since one ionic step corresponds to estep_per_md electronic evolution steps."
                          "\n\n[NOTE] This parameter is only active in RT-TDDFT mode (esolver_type = tddft). It has no effect in ground-state calculations.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.out_freq_td <= 0)
            {
                para.input.out_freq_td = 0; // 0 means no output of info
            }
        };
        read_sync_int(input.out_freq_td);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_elec");
        item.annotation = "print information every few electronic steps";
        item.category = "Output information";
        item.type = "Integer";
        item.description = "Output the charge density (only binary format, controlled by out_chg), wavefunction (controlled by out_wfc_pw) per out_freq_elec electronic iterations. Note that they are always output when converged or reach the maximum iterations scf_nmax.";
        item.default_value = "scf_nmax";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.out_freq_elec <= 0)
            {
                para.input.out_freq_elec = para.input.scf_nmax;
            }
        };
        read_sync_int(input.out_freq_elec);
        this->add_item(item);
    }
    {
        Input_Item item("out_chg");
        item.annotation = "> 0 output charge density for selected electron steps"
                          ", second parameter controls the precision, default is 3.";
        item.category = "Output information";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = R"(The first integer controls whether to output the charge density on real space grids:
* 1: Output the charge density (in Bohr^-3) on real space grids into the density files in the folder OUT.{suffix} too, which can be read in NSCF calculation.

In molecular dynamics simulations, the output frequency is controlled by out_freq_ion.

[NOTE] In the 3.10-LTS version, the file names are SPIN1_CHG.cube and SPIN1_CHG_INI.cube, etc.)";
        item.default_value = "0 3";
        item.unit = "";
        item.availability = "";
			item.read_value = [](const Input_Item& item, Parameter& para) {
				const size_t count = item.get_size();
				if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_chg needs at least 1 value");
				para.input.out_chg[0] = std::stoi(item.str_values[0]);
            para.input.out_chg[1] = 3;
			if (count >= 2) try { para.input.out_chg[1] = std::stoi(item.str_values[1]); }
			catch (const std::invalid_argument&) { /* do nothing */ }
			catch (const std::out_of_range&) {/* do nothing */}
		};
        // reset value in some special case
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
        item.annotation = "output real space potential (with precision 8)";
        item.category = "Output information";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = R"(* 1: Output the total local potential (i.e., local pseudopotential + Hartree potential + XC potential + external electric field (if exists) + dipole correction potential (if exists) + ...) on real space grids (in Ry) into files in the folder OUT.{suffix}. The files are named as:
 * nspin = 1: pots1.cube;
 * nspin = 2: pots1.cube and pots2.cube;
 * nspin = 4: pots1.cube, pots2.cube, pots3.cube, and pots4.cube
* 2: Output the electrostatic potential on real space grids into OUT.{suffix}/pot_es.cube. The Python script named tools/average_pot/aveElecStatPot.py can be used to calculate the average electrostatic potential along the z-axis and outputs it into ElecStaticPot_AVE. Please note that the total local potential refers to the local component of the self-consistent potential, excluding the non-local pseudopotential. The distinction between the local potential and the electrostatic potential is as follows: local potential = electrostatic potential + XC potential.
* 3: Apart from 1, also output the total local potential of the initial charge density. The files are named as:
 * nspin = 1: pots1_ini.cube;
 * nspin = 2: pots1_ini.cube and pots2_ini.cube;
 * nspin = 4: pots1_ini.cube, pots2_ini.cube, pots3_ini.cube, and pots4_ini.cube

The optional second integer controls the output precision. If not provided, the default precision is 8.

In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.

[NOTE] In the 3.10-LTS version, the file names are SPIN1_POT.cube and SPIN1_POT_INI.cube, etc.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
                    const size_t count = item.get_size();
                    if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_pot needs at least 1 value");
                    para.input.out_pot[0] = std::stoi(item.str_values[0]);
                    para.input.out_pot[1] = 8;
                    if (count >= 2) try { para.input.out_pot[1] = std::stoi(item.str_values[1]); }
                    catch (const std::invalid_argument&) { /* do nothing */ }
                    catch (const std::out_of_range&) {/* do nothing */}
            };

        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
                para.input.out_pot[0] = 0;
            }
        };
        sync_intvec(input.out_pot, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_dmk");
        item.annotation = ">0 output density matrix DM(k) for each k-point";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = R"(Whether to output the density matrix for each k-point into files in the folder OUT.${suffix}. The files are named as:
* For gamma only case:
 * nspin = 1 and 4: dm_nao.csr;
 * nspin = 2: dms1_nao.csr and dms2_nao.csr for the two spin channels.
* For multi-k points case:
 * nspin = 1 and 4: dmk1_nao.csr, dmk2_nao.csr, ...;
 * nspin = 2: dmk1s1_nao.csr... and dmk1s2_nao.csr... for the two spin channels.

[NOTE] In the 3.10-LTS version, the parameter is named out_dm and the file names are SPIN1_DM and SPIN2_DM, etc.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
			item.read_value = [](const Input_Item& item, Parameter& para) {
				const size_t count = item.get_size();
				if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_dmk needs at least 1 value");
				para.input.out_dmk[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_dmk[1] = 8;
			if (count >= 2) try { para.input.out_dmk[1] = std::stoi(item.str_values[1]); }
			catch (const std::invalid_argument&) { /* do nothing */ }
			catch (const std::out_of_range&) {/* do nothing */}
			};
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
                para.input.out_dmk[0] = 0;
            }
        };
        sync_intvec(input.out_dmk, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_dmr");
	    item.annotation = "output density matrix DM(R) with respect to lattice vector R (with precision 8)";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = R"(Whether to output the density matrix with Bravias lattice vector R index into files in the folder OUT.${suffix}. The files are named as dmr{s}{spin index}{g}{geometry index}{_nao} + {".csr"}. Here, 's' refers to spin, where s1 means spin up channel while s2 means spin down channel, and the sparse matrix format 'csr' is mentioned in out_mat_hs2. Finally, if out_app_flag is set to false, the file name contains the optional 'g' index for each ionic step that may have different geometries, and if out_app_flag is set to true, the density matrix with respect to Bravias lattice vector R accumulates during ionic steps:
* nspin = 1: dmrs1_nao.csr;
* nspin = 2: dmrs1_nao.csr and dmrs2_nao.csr for the two spin channels.

[NOTE] In the 3.10-LTS version, the parameter is named out_dm1, and the file names are data-DMR-sparse_SPIN0.csr and data-DMR-sparse_SPIN1.csr, etc.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis (multi-k points)";
        item.read_value = [](const Input_Item& item, Parameter& para) {
		    const size_t count = item.get_size();
		    if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_dmr needs at least 1 value");
		    para.input.out_dmr[0] = assume_as_boolean(item.str_values[0]);
		    para.input.out_dmr[1] = 8;
		    if (count >= 2) try { para.input.out_dmr[1] = std::stoi(item.str_values[1]); }
		    catch (const std::invalid_argument&) { /* do nothing */ }
		    catch (const std::out_of_range&) {/* do nothing */}
	    };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.calculation == "get_wf" || para.input.calculation == "get_pchg")
            {
                para.input.out_dmr[0] = 0;
            }
        };
	    item.check_value = [](const Input_Item& item, const Parameter& para) {
		    if (para.sys.gamma_only_local == true && para.input.out_dmr[0])
		    {
			    ModuleBase::WARNING_QUIT("ReadInput", "out_dmr is only valid for multi-k calculation");
		    }
	    };

	    sync_intvec(input.out_dmr, 2, 0);
	    this->add_item(item);
    }
    {
        Input_Item item("out_wfc_pw");
        item.annotation = "output wave functions";
        item.category = "Output information";
        item.type = "Integer";
        item.description = R"(Whether to output the electronic wavefunction coefficients into files and store them in the folder OUT.${suffix}. The files are named as wf{k}{k-point index}{s}{spin index}{g}{geometry index}{e}{electronic iteration index}{_pw} + {".txt"/".dat"}. Here, the s index refers to spin but the label will not show up for non-spin-polarized calculations, where s1 means spin up channel while s2 means spin down channel, and s4 refers to spinor wave functions that contains both spin channels with spin-orbital coupling or noncollinear calculations enabled. For scf or nscf calculations, g index will not appear, but the g index appears for geometry relaxation and molecular dynamics, where one can use the out_freq_ion command to control. To print out the electroinc wave functions every few SCF iterations, use the out_freq_elec command and the e index will appear in the file name.
* 0: no output
* 1: (txt format)
 * non-gamma-only with nspin=1: wfk1_pw.txt, wfk2_pw.txt, ...;
 * non-gamma-only with nspin=2: wfk1s1_pw.txt, wfk1s2_pw.txt, wfk2s1_pw.txt, wfk2s2_pw.txt, ...;
 * non-gamma-only with nspin=4: wfk1s4_pw.txt, wfk2s4_pw.txt, ...;
* 2: (binary format)
 * non-gamma-only with nspin=1: wfk1_pw.dat, wfk2_pw.dat, ...;
 * non-gamma-only with nspin=2: wfk1s1_pw.dat, wfk1s2_pw.dat, wfk2s1_pw.dat, wfk2s2_pw.dat, ...;
 * non-gamma-only with nspin=4: wfk1s4_pw.dat, wfk2s4_pw.dat, ...;

[NOTE] In the 3.10-LTS version, the file names are WAVEFUNC1.dat, WAVEFUNC2.dat, etc.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Output electronic wave functions in plane wave basis, or transform the real-space electronic wave function into plane wave basis (see get_wf option in calculation with NAO basis)";
        read_sync_int(input.out_wfc_pw);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lcao");
        item.annotation = "ouput LCAO wave functions, 0, no output 1: text, 2: binary";
        item.category = "Output information";
        item.type = "Integer";
        item.description = R"(Whether to output the electronic wavefunction coefficients into files and store them in the folder OUT.${suffix}. The files are named as wf{s}{spin index}{k(optional)}{k-point index}{g(optional)}{geometry index1}{_nao} + {".txt"/".dat"}. Here, 's' refers to spin, where s1 means spin up channel while s2 means spin down channel, and 's12' refer to spinor wave functions that contains both spin channels with spin-orbital coupling or noncollinear calculations enabled. In addition, if 'gamma_only' is set to 0, then the optinoal k-point sampling index appears with the k-point index attached to the electronic wave function file names. Finally, if out_app_flag is set to false, the file name contains the optional 'g' index for each ionic step that may have different geometries, and if out_app_flag is set to true, the wave functions accumulate during ionic steps. If the out_app_flag is set to false, a new folder named WFC will be created, and the wave function files will be saved into it.
* 0: no output
* 1: (txt format)
 * gamma-only: wfs1_nao.txt or wfs2_nao.txt, ...;
 * non-gamma-only: wfs1k1_nao.txt or wfs1k2_nao.txt, ...;
* 2: (binary format)
 * gamma-only: wfs1_nao.dat or wfs2_nao.dat, ...;
 * non-gamma-only: wfs1k1_nao.dat or wfs1k2_nao.dat, ....

The corresponding sequence of the orbitals can be seen in Basis Set.

Also controled by out_freq_ion and out_app_flag.

[NOTE] In the 3.10-LTS version, the file names are WFC_NAO_GAMMA1_ION1.txt and WFC_NAO_K1_ION1.txt, etc.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
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
        Input_Item item("out_dos");
        item.annotation = "output energy and dos";
        item.category = "Output information";
        item.type = "Integer";
        item.description = R"(Whether to output the density of states (DOS). For more information, refer to the dos.md.
* 0: no output
* 1: output the density of states (DOS)
 * nspin=1 or 4: doss1g{geom}_{basis}.txt, where geom is the geometry index when cell changes or ions move while basis is either pw or nao.
 * nspin=2: doss1g{geom}_{basis}.txt and doss2g{geom}_{basis}.txt for two spin channles.
* 2: (LCAO) output the density of states (DOS) and the projected density of states (PDOS)
* 3: output the Fermi surface file (fermi.bxsf) in BXSF format that can be visualized by XCrySDen)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("out_ldos");
        item.annotation = "output mode of local density of states, second parameter controls the precision";
        item.category = "Output information";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = R"(Whether to output the local density of states (LDOS), optionally output precision can be set by a second parameter, default is 3.
* 0: no output
* 1: output the partial charge density for given bias (controlled by stm_bias) in cube file format, which can be used to plot scanning tunneling spectroscopys to mimick STM images using the Python script plot.py.
* 2: output LDOS along a line in real space (controlled by ldos_line). Parameters used to control DOS output are also valid for LDOS.
* 3: output both two LDOS modes above.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_ldos should have 1 or 2 values");
            }
            para.input.out_ldos[0] = std::stoi(item.str_values[0]);
            para.input.out_ldos[1] = (count == 2) ? std::stoi(item.str_values[1]) : 3;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_ldos[0] < 0 || para.input.out_ldos[0] > 3)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_ldos should be 0, 1, 2 or 3");
            }
        };
        sync_intvec(input.out_ldos, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_band");
        item.annotation = "output energy and band structure (with precision 8)";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = R"(Whether to output the eigenvalues of the Hamiltonian matrix (in eV) into the running log during electronic iterations and into a file at the end of calculations. The former can be used with the 'out_freq_elec' parameter while the latter option allows the output precision to be set via a second parameter, with a default value of 8. The output file names are:
 * nspin = 1 or 4: eig.txt;
 * nspin = 2: eigs1.txt and eigs2.txt;
 * For more information, refer to the band.md)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("out_proj_band");
        item.annotation = "output projected band structure";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to output the projected band structure. For more information, refer to the band.md";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
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
        Input_Item item("out_stru");
        item.annotation = "output the structure files after each ion step";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to output structure files per ionic step in geometry relaxation calculations into OUT.{istep}_D, where ${istep} is the ionic step.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::vector<std::string> offlist = {"nscf", "get_s", "get_pchg", "get_wf"};
            if (std::find(offlist.begin(), offlist.end(), para.input.calculation) != offlist.end())
            {
                para.input.out_stru = false;
            }
        };
        read_sync_bool(input.out_stru);
        this->add_item(item);
    }
    {
        Input_Item item("out_level");
        item.annotation = "ie(for electrons); i(for ions);";
        item.category = "Output information";
        item.type = "String";
        item.description = R"(Control the output level of information in OUT.{calculation}.log.
* ie: electronic iteration level, which prints useful information for electronic iterations;
* i: geometry relaxation level, which prints some information for geometry relaxations additionally;
* m: molecular dynamics level, which does not print some information for simplicity.)";
        item.default_value = "ie";
        item.unit = "";
        item.availability = "";
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
        add_bool_bcast(sys.out_md_control);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs");
        item.annotation = "output H and S matrix (with precision 8)";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = R"(Whether to print the upper triangular part of the Hamiltonian matrices and overlap matrices for each k-point into files in the directory OUT.${suffix}. The second number controls precision. For more information, please refer to hs_matrix.md. Also controled by out_freq_ion and out_app_flag.
* For gamma only case:
 * nspin = 1: hks1_nao.txt for the Hamiltonian matrix and sks1_nao.txt for the overlap matrix;
 * nspin = 2: hks1_nao.txt and hks2_nao.txt for the Hamiltonian matrix and sks1_nao.txt for the overlap matrix. Note that the code will not output sks2_nao.txt because it is the same as sks1_nao.txt;
 * nspin = 4: hks12_nao.txt for the Hamiltonian matrix and sks12_nao.txt for the overlap matrix.
* For multi-k points case:
 * nspin = 1: hks1k1_nao.txt for the Hamiltonian matrix at the 1st k-point, and sks1k1_nao.txt for the overlap matrix for the 1st k-point, ...;
 * nspin = 2: hks1k1_nao.txt and hks2k1_nao.txt for the two spin channels of the Hamiltonian matrix at the 1st k-point, and sks1k1_nao.txt for the overlap matrix for the 1st k-point. Note that the code will not output sks2k1_nao.txt because it is the same as sks1k1_nao.txt, ...;
 * nspin = 4: hks12k1_nao.txt for the Hamiltonian matrix at the 1st k-point, and sks12k1_nao.txt for the overlap matrix for the 1st k-point, ...;

[NOTE] In the 3.10-LTS version, the file names are data-0-H and data-0-S, etc.)";
        item.default_value = "False 8";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital basis";
			item.read_value = [](const Input_Item& item, Parameter& para) {
				const size_t count = item.get_size();
				if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_mat_hs needs at least 1 value");
				para.input.out_mat_hs[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_mat_hs[1] = 8;
			if (count >= 2) try { para.input.out_mat_hs[1] = std::stoi(item.str_values[1]); }
			catch (const std::invalid_argument&) { /* do nothing */ }
			catch (const std::out_of_range&) {/* do nothing */}
		};
        // reset value in some special case
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
        Input_Item item("out_mat_hs2");
        item.annotation = "output H(R) and S(R) matrix";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print files containing the Hamiltonian matrix and overlap matrix into files in the directory OUT.${suffix}. For more information, please refer to hs_matrix.md."
                          "\n\n[NOTE] In the 3.10-LTS version, the file names are data-HR-sparse_SPIN0.csr and data-SR-sparse_SPIN0.csr, etc.";
        item.default_value = "False";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
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
        Input_Item item("out_mat_tk");
        item.annotation = "output kinetic matrix of electrons T(k)";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = "Whether to print the upper triangular part of the kinetic matrices for each k-point into OUT.${suffix}/tks1ki_nao.txt, where i is the index of k points. One may optionally provide a second parameter to specify the precision."
                          "\n\n[NOTE] In the 3.10-LTS version, the file names are data-TR-sparse_SPIN0.csr, etc.";
        item.default_value = "False [8]";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital basis";
			item.read_value = [](const Input_Item& item, Parameter& para) {
				const size_t count = item.get_size();
				if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_mat_tk needs at least 1 value");
				para.input.out_mat_tk[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_mat_tk[1] = 8;
			if (count >= 2) try { para.input.out_mat_tk[1] = std::stoi(item.str_values[1]); }
			catch (const std::invalid_argument&) { /* do nothing */ }
			catch (const std::out_of_range&) {/* do nothing */}
        };
        sync_intvec(input.out_mat_tk, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_r");
        item.annotation = "output r(R) matrix";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print the matrix representation of the position matrix into a file named rr.csr in the directory OUT.${suffix}. If calculation is set to get_s, the position matrix can be obtained without scf iterations. For more information, please refer to position_matrix.md."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is data-rR-sparse.csr.";
        item.default_value = "False";
        item.unit = "Bohr";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
        read_sync_bool(input.out_mat_r);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.inp.out_mat_r || para.inp.out_mat_hs2 || para.inp.out_mat_t || para.inp.out_mat_dh
                 || para.inp.dm_to_rho)
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
        Input_Item item("out_mat_t");
        item.annotation = "output T(R) matrix";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Generate files containing the kinetic energy matrix. The format will be the same as the Hamiltonian matrix and overlap matrix as mentioned in out_mat_hs2. The name of the files will be trs1_nao.csr and so on. Also controled by out_freq_ion and out_app_flag."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is data-TR-sparse_SPIN0.csr.";
        item.default_value = "False";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
        read_sync_bool(input.out_mat_t);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_dh");
        item.annotation = "output of derivative of H(R) matrix";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print files containing the derivatives of the Hamiltonian matrix. The format will be the same as the Hamiltonian matrix and overlap matrix as mentioned in out_mat_hs2. The name of the files will be dhrxs1_nao.csr, dhrys1_nao.csr, dhrzs1_nao.csr and so on. Also controled by out_freq_ion and out_app_flag."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is data-dHRx-sparse_SPIN0.csr and so on.";
        item.default_value = "False";
        item.unit = "Ry/Bohr";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
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
        Input_Item item("out_mat_ds");
        item.annotation = "output of derivative of S(R) matrix";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print files containing the derivatives of the overlap matrix. The format will be the same as the overlap matrix as mentioned in out_mat_dh. The name of the files will be dsrxs1.csr and so on. Also controled by out_freq_ion and out_app_flag. This feature can be used with calculation get_s."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is data-dSRx-sparse_SPIN0.csr and so on.";
        item.default_value = "False";
        item.unit = "Ry/Bohr";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
        read_sync_bool(input.out_mat_ds);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_mat_ds && para.input.nspin == 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "out_mat_ds is not available for nspin = 4");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_xc");
        item.annotation = "output exchange-correlation matrix in KS-orbital representation";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print the upper triangular part of the exchange-correlation matrices in Kohn-Sham orbital representation: for each k point into files in the directory OUT.i_nao.txt, where {suffix}/vxc_out.dat. If EXX is calculated, the local and EXX part of band energy will also be printed in OUT.{suffix}/vxc_exx_out.dat, respectively. All the vxc_out.dat files contains 3 integers (nk, nspin, nband) followed by nk*nspin*nband lines of energy Hartree and eV."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is k-$k-Vxc and so on.";
        item.default_value = "False";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital (NAO) and NAO-in-PW basis";
        read_sync_bool(input.out_mat_xc);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_xc2");
        item.annotation = "output exchange-correlation matrix in NAO representation";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print the exchange-correlation matrices in numerical orbital representation: in CSR format in the directory OUT.s."
                          "\n\n[NOTE] In the 3.10-LTS version, the file name is Vxc_R_spin$s and so on.";
        item.default_value = "False";
        item.unit = "Ry";
        item.availability = "Numerical atomic orbital (NAO) basis";
        read_sync_bool(input.out_mat_xc2);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_l");
        item.annotation = "output the expectation values of angular momentum operators";
        item.category = "Output information";
        item.type = R"(Boolean \[Integer\](optional))";
        item.description = "Whether to print the expectation value of the angular momentum operator , , and in the basis of the localized atomic orbitals. The files are named OUT.{suffix}_Lx.dat, OUT.{suffix}_Ly.dat, and OUT.{suffix}_Lz.dat. The second integer controls the precision of the output.";
        item.default_value = "False 8";
        item.unit = "";
        item.availability = "Numerical atomic orbital (NAO) basis";
			item.read_value = [](const Input_Item& item, Parameter& para) {
				const size_t count = item.get_size();
				if (count < 1) ModuleBase::WARNING_QUIT("ReadInput", "out_mat_l needs at least 1 value");
				para.input.out_mat_l[0] = assume_as_boolean(item.str_values[0]);
            para.input.out_mat_l[1] = 8;
			if (count >= 2) try { para.input.out_mat_l[1] = std::stoi(item.str_values[1]); }
			catch (const std::invalid_argument&) { /* do nothing */ }
			catch (const std::out_of_range&) {/* do nothing */}
        };
        sync_intvec(input.out_mat_l, 2, 0);
        this->add_item(item);
    }
    {
        Input_Item item("out_xc_r");
        item.annotation = "if >=0, output the derivatives of exchange correlation in realspace, second parameter controls the precision";
        item.category = "Output information";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = R"(The first integer controls whether to output the exchange-correlation (in Bohr^-3) on real space grids using Libxc to folder OUT.${suffix}:
* 0: rho, amag, sigma, exc
* 1: vrho, vsigma
* 2: v2rho2, v2rhosigma, v2sigma2
* 3: v3rho3, v3rho2sigma, v3rhosigma2, v3sigma3
* 4: v4rho4, v4rho3sigma, v4rho2sigma2, v4rhosigma3, v4sigma4 The meaning of the files is presented in Libxc

The second integer controls the precision of the charge density output, if not given, will use 3 as default.

The circle order of the charge density on real space grids is: x is the outer loop, then y and finally z (z is moving fastest).)";
        item.default_value = "-1 3";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count >= 1)
            {
                para.input.out_xc_r[0] = std::stoi(item.str_values[0]);
            }
            if (count >= 2)
            {
                para.input.out_xc_r[1] = std::stoi(item.str_values[1]);
            }
        };
        // check value
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.out_xc_r[0] >= 0)
            {
#ifndef USE_LIBXC
                ModuleBase::WARNING_QUIT("ReadInput", "INPUT out_xc_r is only aviailable with Libxc");
#endif
            }
        };
        sync_intvec(input.out_xc_r, 2, -1);
        this->add_item(item);
    }
    {
        Input_Item item("out_eband_terms");
        item.annotation = "output the band energy terms separately";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print the band energy terms separately in the file OUT.{term}_out.dat. The terms include the kinetic, pseudopotential (local + nonlocal), Hartree and exchange-correlation (including exact exchange if calculated).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_bool(input.out_eband_terms);
        this->add_item(item);
    }
    {
        Input_Item item("out_mul");
        item.annotation = "mulliken charge or not";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print the Mulliken population analysis result into OUT.${suffix}/mulliken.txt. In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
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
        Input_Item item("out_app_flag");
        item.annotation = "whether output r(R), H(R), S(R), T(R), and dH(R) "
                          "matrices in an append manner during MD";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to output r(R), H(R), S(R), T(R), dH(R), dS(R), and wfc matrices in an append manner during molecular dynamics calculations. Check input parameters out_mat_r, out_mat_hs2, out_mat_t, out_mat_dh, out_mat_hs and out_wfc_lcao for more information.";
        item.default_value = "true";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis (not gamma-only algorithm)";
        read_sync_bool(input.out_app_flag);
        this->add_item(item);
    }
    {
        Input_Item item("out_ndigits");
        item.annotation = "the length of decimal part of output data";
        item.category = "Output information";
        item.type = "Integer";
        item.description = "Controls the length of decimal part of output data, such as charge density, Hamiltonian matrix, Overlap matrix and so on.";
        item.default_value = "8";
        item.unit = "";
        item.availability = "out_mat_hs 1 case presently.";
        read_sync_int(input.out_ndigits);
        this->add_item(item);
    }
    {
        Input_Item item("out_element_info");
        item.annotation = "output (projected) wavefunction of each element";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Whether to print element information into files in the directory OUT.{element_label}, including pseudopotential and orbital information of the element (in atomic Ryberg units).";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_element_info);
        this->add_item(item);
    }
    {
        Input_Item item("restart_save");
        item.annotation = "print to disk every step for restart";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = R"(Whether to save charge density files per ionic step, which are used to restart calculations. According to the value of read_file_dir:
* auto: These files are saved in folder OUT.{read_file_dir}/restart/.

If EXX(exact exchange) is calculated (i.e. dft_fuctional==hse/hf/pbe0/scan0 or rpa==True), the Hexx(R) files for each processor will also be saved in the above folder, which can be read in EXX calculation with restart_load==True.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_bool(input.restart_save);
        this->add_item(item);
    }
    {
        Input_Item item("rpa");
        item.annotation = "true:generate output files used in rpa calculation; "
                          "false:(default)";
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Generate output files used in rpa calculations."
                          "\n\n[NOTE] If symmetry is set to 1, additional files containing the necessary information for "
                          "exploiting symmetry in the subsequent rpa calculation will be output: "
                          "irreducible_sector.txt, symrot_k.txt and symrot_R.txt.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.rpa);
        this->add_item(item);
    }
    {
        Input_Item item("out_pchg");
        item.annotation = "specify the bands to be calculated for the partial (band-decomposed) charge densities";
        item.category = "Output information";
        item.type = "String";
        item.description = R"(Specifies the electronic states to calculate the charge densities with state index for, using a space-separated string of 0s and 1s. Each digit in the string corresponds to a state, starting from the first state. A 1 indicates that the charge density should be calculated for that state, while a 0 means the state will be ignored. The parameter allows a compact and flexible notation (similar to ocp_set), for example the syntax 1 4*0 5*1 0 is used to denote the selection of states: 1 means calculate for the first state, 4*0 skips the next four states, 5*1 means calculate for the following five states, and the final 0 skips the next state. It's essential that the total count of states does not exceed the total number of states (nbands); otherwise, it results in an error, and the process exits. The input string must contain only numbers and the asterisk (*) for repetition, ensuring correct format and intention of state selection. The outputs comprise multiple .cube files following the naming convention pchgi[state]s[spin]k[kpoint].cube.)";
        item.default_value = "none";
        item.unit = "";
        item.availability = "For both PW and LCAO. When basis_type = lcao, used when calculation = get_pchg.";
        item.read_value
            = [](const Input_Item& item, Parameter& para) { parse_expression(item.str_values, para.input.out_pchg); };
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
        item.category = "Output information";
        item.type = "String";
        item.description = "Specifies the electronic states to calculate the real-space wave function modulus (norm, or known as the envelope function) with state index. The syntax and state selection rules are identical to out_pchg, but the output is the norm of the wave function. The outputs comprise multiple .cube files following the naming convention wfi[state]s[spin]k[kpoint].cube.";
        item.default_value = "none";
        item.unit = "";
        item.availability = "For both PW and LCAO. When basis_type = lcao, used when calculation = get_wf.";
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
        item.category = "Output information";
        item.type = "String";
        item.description = "Specifies the electronic states to calculate the real and imaginary parts of the wave function with state index. The syntax and state selection rules are identical to out_pchg, but the output contains both the real and imaginary components of the wave function. The outputs comprise multiple .cube files following the naming convention wfi[state]s[spin]k[kpoint][re/im].cube.";
        item.default_value = "none";
        item.unit = "";
        item.availability = "For both PW and LCAO. When basis_type = lcao, used when calculation = get_wf.";
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
        item.category = "Output information";
        item.type = "Boolean";
        item.description = "Specifies whether to write the partial charge densities for all k-points to individual files or merge them. Warning: Enabling symmetry may produce unwanted results due to reduced k-point weights and symmetry operations in real space. Therefore when calculating partial charge densities, if you are not sure what you want exactly, it is strongly recommended to set symmetry = -1. It is noteworthy that your symmetry setting should remain the same as that in the SCF procedure.";
        item.default_value = "false";
        item.unit = "";
        item.availability = "For both PW and LCAO. When basis_type = pw, used if out_pchg is set. When basis_type = lcao, used only when calculation = get_pchg and gamma_only = 0.";
        read_sync_bool(input.if_separate_k);
        this->add_item(item);
    }
    {
        Input_Item item("out_elf");
        item.annotation = "> 0 output electron localization function (ELF) for selected electron steps"
                          ", second parameter controls the precision, default is 3.";
        item.category = "Output information";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = R"(Whether to output the electron localization function (ELF) in the folder `OUT.${suffix}`. The files are named as
* nspin = 1:
    * elf.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i}|^2} - \frac{|\nabla\rho|^2}{8\rho}}{\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}}$;
* nspin = 2:
    * elf1.cube, elf2.cube: ${\rm{ELF}}_\sigma = \frac{1}{1+\chi_\sigma^2}$, $\chi_\sigma = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i,\sigma}|^2} - \frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}$;
    * elf.cube: ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i,\sigma}{f_i |\nabla\psi_{i,\sigma}|^2} - \sum_{\sigma}{\frac{|\nabla\rho_\sigma|^2}{8\rho_\sigma}}}{\sum_{\sigma}{\frac{3}{10}(6\pi^2)^{2/3}\rho_\sigma^{5/3}}}$;
* nspin = 4 (noncollinear):
    * elf.cube: ELF for total charge density, ${\rm{ELF}} = \frac{1}{1+\chi^2}$, $\chi = \frac{\frac{1}{2}\sum_{i}{f_i |\nabla\psi_{i}|^2} - \frac{|\nabla\rho|^2}{8\rho}}{\frac{3}{10}(3\pi^2)^{2/3}\rho^{5/3}}$

The second integer controls the precision of the kinetic energy density output, if not given, will use 3 as default. For purpose restarting from this file and other high-precision involved calculation, recommend to use 10.

In molecular dynamics calculations, the output frequency is controlled by out_freq_ion.)";
        item.default_value = "0 3";
        item.unit = "";
        item.availability = "Only for Kohn-Sham DFT and Orbital Free DFT.";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count >= 1)
            {
                para.input.out_elf[0] = std::stoi(item.str_values[0]);
            }
            if (count >= 2)
            {
                para.input.out_elf[1] = std::stoi(item.str_values[1]);
            }
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
    {
        // refactored from the removal of wannier input file, ISSUE 6469
        Input_Item item("out_spillage");
        item.annotation = "output spillage of wavefunctions. This parameter only accepts 0 or 2.";
        item.category = "Output information";
        item.type = "Integer";
        item.description = "This output is only intentively needed by the ABACUS numerical atomic orbital generation workflow. This parameter is used to control whether to output the overlap integrals between truncated spherical Bessel functions (TSBFs) and plane-wave basis expanded wavefunctions (named as OVERLAP_Q), and between TSBFs (named as OVERLAP_Sq), also their first order derivatives. The output files are named starting with orb_matrix. A value of 2 would enable the output.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Only for Kohn-Sham DFT with plane-wave basis.";
        read_sync_int(input.out_spillage);
        this->add_item(item);
    }
    {
        Input_Item item("out_dipole");
        item.annotation = "output dipole or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(* True: Output electric dipole moment.
* False: Do not output electric dipole moment.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_dipole);
        this->add_item(item);
    }
    {
        Input_Item item("out_current");
        item.annotation = "output current or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(* 0: Do not output current.
* 1: Output current using the two-center integral, faster.
* 2: Output current using the matrix commutation, more precise.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.out_current);
        this->add_item(item);
    }
    {
        Input_Item item("out_current_k");
        item.annotation = "output current for each k";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(* True: Output current for each k-points separately.
* False: Output current in total.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_current_k);
        this->add_item(item);
    }
    {
        Input_Item item("out_efield");
        item.annotation = "output dipole or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Whether to output the electric field data to files. When enabled, writes real-time electric field values (unit: V/A) into files named efield_[num].txt, where [num] is the sequential index of the electric field ranges from 0 to N-1 for N configured fields. It is noteworthy that the field type sequence follows td_ttype, while the direction sequence follows td_vext_dire.
* True: Output electric field.
* False: Do not output electric field.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_efield);
        this->add_item(item);
    }
    {
        Input_Item item("out_vecpot");
        item.annotation = "output TDDFT vector potential or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Output vector potential or not (unit: a.u.).
* True: Output vector potential into file At.dat.
* False: Do not output vector potential.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_vecpot);
        this->add_item(item);
    }
    {
        // recover the functionality of test_symmetry by introducing a new keyword "out_symm_mat"
        // the "out_symm_mat" keyword will be a
        Input_Item item("cal_symm_repr");
        item.annotation = "output matrix representation of symmetry operation into running log file"
                          " > 0 output the matrix representation of symmetry operation "
                          ", the second parameter controls the precision, default is 3.";
        item.category = "System variables";
        item.type = R"(Integer \[Integer\](optional))";
        item.description = "Whether to print the matrix representation of symmetry operation to running log file. If the first value is given as 1, then all matrix representations will be printed. The second optional parameter controls the precision (number of digits) to print, default is 3, which is enough for a quick check.";
        item.default_value = "1 3";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count >= 1)
            {
                para.input.cal_symm_repr[0] = std::stoi(item.str_values[0]);
            }
            if (count >= 2)
            {
                para.input.cal_symm_repr[1] = std::stoi(item.str_values[1]);
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.cal_symm_repr[0] < 0 || para.input.cal_symm_repr[0] > 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "cal_symm_repr should be 0 or 1");
            }
        };
        sync_intvec(input.cal_symm_repr, 2, 0);
        this->add_item(item);
    }
    {
        // refactored from the removal of wannier input file, ISSUE 6469
        Input_Item item("spillage_outdir");
        item.annotation = "output directory for spillage of wavefunctions.";
        item.category = "Input files";
        item.type = "String";
        item.description = "The directory to save the spillage files.";
        item.default_value = "\"./\"";
        item.unit = "";
        item.availability = "Used only for plane wave basis set.";
        read_sync_string(input.spillage_outdir);
        this->add_item(item);
    }
}
} // namespace ModuleIO
