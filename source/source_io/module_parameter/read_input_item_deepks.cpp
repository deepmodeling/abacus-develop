#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_io/module_parameter/parameter.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_deepks()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("deepks_out_labels");
        item.annotation = ">0 compute descriptor for deepks. 1 used during training, 2 used for label production";
        item.category = "DeePKS";
        item.type = "Integer";
        item.description = R"(Print labels and descriptors for DeePKS in OUT.${suffix}. The names of these files start with "deepks".
* 0 : No output.
* 1 : Output intermediate files needed during DeePKS training.
* 2 : Output target labels for label preperation. The label files are named as deepks_<property>.npy or deepks_<property>.csr, where the units and formats are the same as label files <property>.npy or <property>.csr required for training, except that the first dimension (nframes) is excluded. System structrue files are also given in deepks_atom.npy and deepks_box.npy in the unit of Bohr, which means lattice_constant should be set to 1 when training.

[NOTE] When deepks_out_labels equals 1, the path of a numerical descriptor (an orb file) is needed to be specified under the NUMERICAL_DESCRIPTOR tag in the STRU file. This is not needed when deepks_out_labels equals 2.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_int(input.deepks_out_labels);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_freq_elec");
        item.annotation = ">0 frequency of electronic iteration to output descriptors and labels for deepks.";
        item.category = "DeePKS";
        item.type = "Integer";
        item.description = "When deepks_out_freq_elec is greater than 0, print labels and descriptors for DeePKS in OUT.${suffix}/DeePKS_Labels_Elec per deepks_out_freq_elec electronic iterations, with suffix _e* to distinguish different steps. Often used with deepks_out_labels equals 1.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_int(input.deepks_out_freq_elec);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.deepks_out_freq_elec < 0)
            {
                para.input.deepks_out_freq_elec = 0;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_out_freq_elec > 0 && para.input.deepks_out_base == "none")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "to use deepks_out_freq_elec, please set deepks_out_base ");
            }            
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_base");
        item.annotation = "base functional for output files, with dft_functional as target functional";
        item.category = "DeePKS";
        item.type = "String";
        item.description = "Print labels and descriptors calculated by base functional ( determined by deepks_out_base ) and target functional ( determined by dft_functional ) for DeePKS in per deepks_out_freq_elec electronic iterations. The SCF process, labels and descriptors output of the target functional are all consistent with those when the target functional is used alone. The only additional output under this configuration is the labels of the base functional. Often used with deepks_out_labels equals 1.";
        item.default_value = "None";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis and deepks_out_freq_elec is greater than 0";
        read_sync_string(input.deepks_out_base);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_out_base != "none" && para.input.deepks_out_labels == 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "to use deepks_out_base, please set deepks_out_labels > 0 ");
            }
            if (para.input.deepks_out_base != "none" && para.input.deepks_bandgap > 0 )
            {
                ModuleBase::WARNING_QUIT("ReadInput", "outputting bandgap labels during electronic steps is not implemented yet ");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_scf");
        item.annotation = ">0 add V_delta to Hamiltonian";
        item.category = "DeePKS";
        item.type = "Boolean";
        item.description = "perform self-consistent field iteration in DeePKS method"
                          "\n\n[NOTE] A trained, traced model file is needed.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_bool(input.deepks_scf);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
#ifndef __MLALGO
            if (para.input.deepks_scf || para.input.deepks_out_labels || para.input.deepks_bandgap
                || para.input.deepks_v_delta)
            {
                ModuleBase::WARNING_QUIT("Input_conv", "please compile with DeePKS");
            }
#endif
            // if (!para.input.deepks_scf && para.input.deepks_out_labels == 1)
            // {
            //     ModuleBase::WARNING_QUIT("Input_conv", "deepks_out_labels = 1 requires deepks_scf = 1");
            // }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_equiv");
        item.annotation = "whether to use equivariant version of DeePKS";
        item.category = "DeePKS";
        item.type = "Boolean";
        item.description = "whether to use equivariant version of DeePKS"
                          "\n\n[NOTE] The equivariant version of DeePKS-kit is still under development, "
                          "so this feature is currently only intended for internal usage.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_bool(input.deepks_equiv);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.deepks_equiv && para.input.deepks_bandgap)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "equivariant version of DeePKS is not implemented yet");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_model");
        item.annotation = "file dir of traced pytorch model: 'model.ptg";
        item.category = "DeePKS";
        item.type = "String";
        item.description = "the path of the trained, traced neural network model file generated by deepks-kit";
        item.default_value = "None";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis and deepks_scf is true";
        read_sync_string(input.deepks_model);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_lmax");
        item.annotation = "lmax used in generating spherical bessel functions";
        item.category = "DeePKS";
        item.type = "Integer";
        item.description = "the maximum angular momentum of the Bessel functions generated as the projectors in DeePKS - NOte: To generate such projectors, set calculation type to gen_bessel in ABACUS. See also calculation.";
        item.default_value = "2";
        item.unit = "";
        item.availability = "gen_bessel calculation";
        read_sync_int(input.bessel_descriptor_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
        item.category = "DeePKS";
        item.type = "String";
        item.description = "energy cutoff of Bessel functions";
        item.default_value = "same as ecutwfc";
        item.unit = "Ry";
        item.availability = "gen_bessel calculation";
        read_sync_string(input.bessel_descriptor_ecut);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.bessel_descriptor_ecut == "default")
            {
                para.input.bessel_descriptor_ecut = std::to_string(para.input.ecutwfc);
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (std::stod(para.input.bessel_descriptor_ecut) < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_descriptor_ecut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_tolerence");
        item.annotation = "tolerence for spherical bessel root";
        item.category = "DeePKS";
        item.type = "Real";
        item.description = "tolerance for searching the zeros of Bessel functions";
        item.default_value = "1.0e-12";
        item.unit = "";
        item.availability = "gen_bessel calculation";
        read_sync_double(input.bessel_descriptor_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
        item.category = "DeePKS";
        item.type = "Real";
        item.description = "cutoff radius of Bessel functions";
        item.default_value = "6.0";
        item.unit = "Bohr";
        item.availability = "gen_bessel calculation";
        read_sync_double(input.bessel_descriptor_rcut);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.bessel_descriptor_rcut < 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "bessel_descriptor_rcut must >= 0");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_smooth");
        item.annotation = "spherical bessel smooth or not";
        item.category = "DeePKS";
        item.type = "Boolean";
        item.description = "smooth the Bessel functions at radius cutoff";
        item.default_value = "False";
        item.unit = "";
        item.availability = "gen_bessel calculation";
        read_sync_bool(input.bessel_descriptor_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_sigma");
        item.annotation = "sphereical bessel smearing_sigma";
        item.category = "DeePKS";
        item.type = "Real";
        item.description = "smooth parameter at the cutoff radius of projectors";
        item.default_value = "0.1";
        item.unit = "Bohr";
        item.availability = "gen_bessel calculation";
        read_sync_double(input.bessel_descriptor_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_bandgap");
        item.annotation = ">0 for bandgap label";
        item.category = "DeePKS";
        item.type = "Integer";
        item.description = R"(include bandgap label for DeePKS training
* 0: Don't include bandgap label
* 1: Include target bandgap label (see deepks_band_range for more details)
* 2: Include multiple bandgap label (see deepks_band_range for more details)
* 3: Used for systems containing H atoms. Here HOMO is defined as the max occupation except H atoms and the bandgap label is the energy between HOMO and (HOMO + 1))";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis and deepks_scf is true";
        read_sync_int(input.deepks_bandgap);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_bandgap < 0 || para.input.deepks_bandgap > 3)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "deepks_bandgap must be integer in [0, 3]");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_band_range");
        item.annotation = "(int, int) range of bands for bandgap label";
        item.category = "DeePKS";
        item.type = "Int*2";
        item.description = R"(The first value should not be larger than the second one and the meaning differs in different cases below
* deepks_bandgap is 1: Bandgap label is the energy between LUMO + deepks_band_range[0] and LUMO + deepks_band_range[1]. If not set, it will calculate energy between HOMO and LUMO states.
* deepks_bandgap is 2: Bandgap labels are energies between HOMO and all states in range [LUMO + deepks_band_range[0], LUMO + deepks_band_range[1]] (Thus there are deepks_band_range[1] - deepks_band_range[0] + 1 bandgaps in total). If HOMO is included in the setting range, it will be ignored since it will always be zero and has no valuable messages (deepks_band_range[1] - deepks_band_range[0] bandgaps in this case). NOTICE: The set range can be greater than, less than, or include the value of HOMO. In the bandgap label, we always calculate the energy of the state in the set range minus the energy of HOMO state, so the bandgap can be negative if the state is lower than HOMO.)";
        item.default_value = "-1 0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis, deepks_scf is true, and deepks_bandgap is 1 or 2";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.deepks_band_range[0] = std::stod(item.str_values[0]);
            para.input.deepks_band_range[1] = std::stod(item.str_values[1]);
        };
        sync_intvec(input.deepks_band_range, 2, 0);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_bandgap == 1)
            {
                if (para.input.deepks_band_range[0] >= para.input.deepks_band_range[1])
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "deepks_band_range[0] must be smaller than deepks_band_range[1] for deepks_bandgap = 1.");
                }
            }
            else if (para.input.deepks_bandgap == 2)
            {
                if (para.input.deepks_band_range[0] > para.input.deepks_band_range[1])
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "deepks_band_range[0] must be no more than deepks_band_range[1] for deepks_bandgap = 2.");
                }
            }
            else
            {
                if (para.input.deepks_band_range[0] != -1 || para.input.deepks_band_range[1] != 0)
                {
                    ModuleBase::WARNING("ReadInput", "deepks_band_range is used for deepks_bandgap = 1/2. Ignore its setting for other cases.");
                }
            } 
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_v_delta");
        item.annotation = "!=0 for v_delta/v_delta_R label. when output, 1 for vdpre, 2 for phialpha and grad_evdm (can save memory )"
                          " -1 for vdrpre, -2 for phialpha_r and gevdm (can save memory)";
        item.category = "DeePKS";
        item.type = "Integer";
        item.description = R"(Include V_delta/V_delta_R (Hamiltonian in k/real space) label for DeePKS training. When deepks_out_labels is true and deepks_v_delta > 0 (k space), ABACUS will output deepks_hbase.npy, deepks_vdelta.npy and deepks_htot.npy(htot=hbase+vdelta). When deepks_out_labels is true and deepks_v_delta < 0 (real space), ABACUS will output deepks_hrtot.csr, deepks_hrdelta.csr. Some more files output for different settings. NOTICE: To match the unit Normally used in DeePKS, the unit of Hamiltonian in k space is Hartree. However, currently in R space the unit is still Ry.
* deepks_v_delta = 1: deepks_vdpre.npy, which is used to calculate V_delta during DeePKS training.
* deepks_v_delta = 2: deepks_phialpha.npy and deepks_gevdm.npy, which can be used to calculate deepks_vdpre.npy. A recommanded method for memory saving.
* deepks_v_delta = -1: deepks_vdrpre.npy, which is used to calculate V_delta_R during DeePKS training.
* deepks_v_delta = -2: deepks_phialpha_r.npy and deepks_gevdm.npy, which can be used to calculate deepks_vdrpre.npy. A recommanded method for memory saving.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "Numerical atomic orbital basis";
        read_sync_int(input.deepks_v_delta);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_v_delta < -2 || para.input.deepks_v_delta > 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "deepks_v_delta must be integer in [-2, 2]");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_unittest");
        item.annotation = "if set 1, prints intermediate quantities that shall "
                          "be used for making unit test";
        item.category = "DeePKS";
        item.type = "Boolean";
        item.description = "generate files for constructing DeePKS unit test"
                          "\n\n[NOTE] Not relevant when running actual calculations. "
                          "When set to 1, ABACUS needs to be run with only 1 process.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.deepks_out_unittest);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.deepks_out_unittest)
            {
                para.input.deepks_out_labels = 1;
                para.input.deepks_scf = true;
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_out_unittest)
            {
                if (para.input.cal_force != 1 || para.input.cal_stress != 1)
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "force and stress are required in generating deepks unittest");
                }
            }
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO
