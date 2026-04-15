#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_postprocess()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("dos_edelta_ev");
        item.annotation = "delta energy for dos";
        item.category = "Density of states";
        item.type = "Real";
        item.description = "The step size in writing Density of States (DOS)";
        item.default_value = "0.01";
        item.unit = "eV";
        item.availability = "";
        read_sync_double(input.dos_edelta_ev);
        this->add_item(item);
    }
    {
        Input_Item item("dos_sigma");
        item.annotation = "gauss b coefficeinet(default=0.07)";
        item.category = "Density of states";
        item.type = "Real";
        item.description = "The width of the Gaussian factor when obtaining smeared Density of States (DOS)";
        item.default_value = "0.07";
        item.unit = "eV";
        item.availability = "";
        read_sync_double(input.dos_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("dos_scale");
        item.annotation = "scale dos range by";
        item.category = "Density of states";
        item.type = "Real";
        item.description = "Defines the energy range of DOS output as (emax-emin)*(1+dos_scale), centered at (emax+emin)/2. This parameter will be used when dos_emin and dos_emax are not set.";
        item.default_value = "0.01";
        item.unit = "eV";
        item.availability = "";
        read_sync_double(input.dos_scale);
        this->add_item(item);
    }
    // DOS
    {
        Input_Item item("dos_emin_ev");
        item.annotation = "minimal range for dos";
        item.category = "Density of states";
        item.type = "Real";
        item.description = R"(The minimal range for Density of States (DOS)
* If set, "dos_scale" will be ignored.)";
        item.default_value = "Minimal eigenenergy of";
        item.unit = "eV";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emin_ev = doublevalue;
            para.sys.dos_setemin = true;
        };
        sync_double(input.dos_emin_ev);
        add_bool_bcast(sys.dos_setemin);
        this->add_item(item);
    }
    {
        Input_Item item("dos_emax_ev");
        item.annotation = "maximal range for dos";
        item.category = "Density of states";
        item.type = "Real";
        item.description = R"(The maximal range for Density of States (DOS)
* If set, "dos_scale" will be ignored.)";
        item.default_value = "Maximal eigenenergy of";
        item.unit = "eV";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emax_ev = doublevalue;
            para.sys.dos_setemax = true;
        };
        sync_double(input.dos_emax_ev);
        add_bool_bcast(sys.dos_setemax);
        this->add_item(item);
    }
    {
        Input_Item item("dos_nche");
        item.annotation = "orders of Chebyshev expansions for dos";
        item.category = "Density of states";
        item.type = "Integer";
        item.description = "The order of Chebyshev expansions when using Stochastic Density Functional Theory (SDFT) to calculate DOS.";
        item.default_value = "100";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.dos_nche);
        this->add_item(item);
    }
    {
        Input_Item item("stm_bias");
        item.annotation = "bias voltage used to calculate ldos";
        item.category = "Density of states";
        item.type = "Real Real(optional) Integer(optional)";
        item.description = R"(The bias voltage used to calculate local density of states to simulate scanning tunneling microscope, see details in out_ldos. When using three parameters:

* The first parameter specifies the initial bias voltage value.
* The second parameter defines the voltage increment (step size between consecutive bias values).
* The third parameter determines the total number of voltage points)";
        item.default_value = "1.0";
        item.unit = "V";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 1 && count != 3)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "stm_bias should have 1 or 3 values");
            }
            para.input.stm_bias[0] = std::stod(item.str_values[0]);
            para.input.stm_bias[1] = (count == 3) ? std::stod(item.str_values[1]) : 0.1;
            para.input.stm_bias[2] = (count == 3) ? std::stod(item.str_values[2]) : 1;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.stm_bias[2] <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "stm_bias[2] should be greater than 0");
            }
            if (para.input.stm_bias[1] == 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "stm_bias[1] should be nonzero");
            }
        };
        sync_doublevec(input.stm_bias, 3, 0);
        this->add_item(item);
    }
    {
        Input_Item item("ldos_line");
        item.annotation = "start and end point of the line (direct coordinates) and number of points";
        item.category = "Density of states";
        item.type = "Real*6 Integer(optional)";
        item.description = "Specify the path of the three-dimensional space and display LDOS in the form of a two-dimensional color chart, see details in out_ldos. The first three paramenters are the direct coordinates of the start point, the next three paramenters are the direct coordinates of the end point, and the final one is the number of points along the path, whose default is 100.";
        item.default_value = "0.0 0.0 0.0 0.0 0.0 1.0 100";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            const size_t count = item.get_size();
            if (count != 6 && count != 7)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ldos_line should have 6 or 7 values");
            }
            for (int i = 0; i < 6; ++i)
            {
                para.input.ldos_line[i] = std::stod(item.str_values[i]);
            }
            para.input.ldos_line[6] = (count == 7) ? std::stoi(item.str_values[6]) : 100;
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.ldos_line[6] <= 0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "ldos_line[6] should be greater than 0");
            }
        };
        sync_doublevec(input.ldos_line, 7, 0);
        this->add_item(item);
    }

    // Electronic Conductivity
    {
        Input_Item item("cal_cond");
        item.annotation = "calculate electronic conductivities";
        item.category = "Electronic conductivities";
        item.type = "Boolean";
        item.description = "Whether to calculate electronic conductivities.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "basis_type = pw";
        read_sync_bool(input.cal_cond);
        this->add_item(item);
    }
    {
        Input_Item item("cond_che_thr");
        item.annotation = "control the error of Chebyshev expansions for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Real";
        item.description = "Control the error of Chebyshev expansions for conductivities.";
        item.default_value = "1e-8";
        item.unit = "";
        item.availability = "esolver_type = sdft";
        read_sync_double(input.cond_che_thr);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dw");
        item.annotation = "frequency interval for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Real";
        item.description = "Frequency interval () for frequency-dependent conductivities.";
        item.default_value = "0.1";
        item.unit = "eV";
        item.availability = "basis_type = pw";
        read_sync_double(input.cond_dw);
        this->add_item(item);
    }
    {
        Input_Item item("cond_wcut");
        item.annotation = "cutoff frequency (omega) for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Real";
        item.description = "Cutoff frequency for frequency-dependent conductivities.";
        item.default_value = "10.0";
        item.unit = "eV";
        item.availability = "basis_type = pw";
        read_sync_double(input.cond_wcut);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dt");
        item.annotation = "t interval to integrate Onsager coefficiencies";
        item.category = "Electronic conductivities";
        item.type = "Real";
        item.description = "Time interval () to integrate Onsager coefficients.";
        item.default_value = "0.02";
        item.unit = "a.u.";
        item.availability = "basis_type = pw";
        read_sync_double(input.cond_dt);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dtbatch");
        item.annotation = "exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion";
        item.category = "Electronic conductivities";
        item.type = "Integer";
        item.description = R"(exp(iH\dt\cond_dtbatch) is expanded with Chebyshev expansion to calculate conductivities. It is faster but costs more memory.
* If cond_dtbatch = 0: Autoset this parameter to make expansion orders larger than 100.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "esolver_type = sdft";
        read_sync_int(input.cond_dtbatch);
        this->add_item(item);
    }
    {
        Input_Item item("cond_smear");
        item.annotation = "Smearing method for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Integer";
        item.description = R"(Smearing method for conductivities
* 1: Gaussian smearing
* 2: Lorentzian smearing)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.cond_smear);
        this->add_item(item);
    }
    {
        Input_Item item("cond_fwhm");
        item.annotation = "FWHM for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Real";
        item.description = "FWHM for conductivities. For Gaussian smearing, ; for Lorentzian smearing, .";
        item.default_value = "0.4";
        item.unit = "eV";
        item.availability = "basis_type = pw";
        read_sync_double(input.cond_fwhm);
        this->add_item(item);
    }
    {
        Input_Item item("cond_nonlocal");
        item.annotation = "Nonlocal effects for conductivities";
        item.category = "Electronic conductivities";
        item.type = "Boolean";
        item.description = R"(Whether to consider nonlocal potential correction when calculating velocity matrix .
* True: .
* False: .)";
        item.default_value = "True";
        item.unit = "";
        item.availability = "basis_type = pw";
        read_sync_bool(input.cond_nonlocal);
        this->add_item(item);
    }

    // berry_wannier
    {
        Input_Item item("berry_phase");
        item.annotation = "calculate berry phase or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Controls the calculation of Berry phase
* true: Calculate Berry phase.
* false: Do not calculate Berry phase.)";
        item.default_value = "false";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.berry_phase);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.berry_phase)
            {
                if (para.input.basis_type != "pw" && para.input.basis_type != "lcao")
                {
                    ModuleBase::WARNING_QUIT("ReadInput",
                                             "calculate berry phase, please "
                                             "set basis_type = pw or lcao");
                }
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set calculation = nscf");
                }
                if (!(para.input.gdir == 1 || para.input.gdir == 2 || para.input.gdir == 3))
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set gdir = 1 or 2 or 3");
                }
                if (para.input.symmetry != "-1")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "calculate berry phase, please set symmetry = -1");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("gdir");
        item.annotation = "calculate the polarization in the direction of the "
                          "lattice vector";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Integer";
        item.description = R"(The direction of the polarization in the lattice vector for Berry phase calculation
* 1: Calculate the polarization in the direction of the lattice vector a_1 defined in the STRU file.
* 2: Calculate the polarization in the direction of the lattice vector a_2 defined in the STRU file.
* 3: Calculate the polarization in the direction of the lattice vector a_3 defined in the STRU file.)";
        item.default_value = "3";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.gdir);
        this->add_item(item);
    }
    {
        Input_Item item("towannier90");
        item.annotation = "use wannier90 code interface or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Controls the generation of files for the Wannier90 code.
* 1: Generate files for the Wannier90 code.
* 0: Do not generate files for the Wannier90 code.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.towannier90);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.towannier90)
            {
                if (para.input.calculation != "nscf")
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "to use towannier90, please set calculation = nscf");
                }
                if (para.input.nspin == 2)
                {
                    if (para.input.wannier_spin != "up" && para.input.wannier_spin != "down")
                    {
                        ModuleBase::WARNING_QUIT("ReadInput",
                                                 "to use towannier90, please set wannier_spin = up "
                                                 "or down");
                    }
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("nnkpfile");
        item.annotation = "the wannier90 code nnkp file name";
        item.category = "Berry phase and wannier90 interface";
        item.type = "String";
        item.description = "The file name generated when running \"wannier90 -pp ...\" command";
        item.default_value = "seedname.nnkp";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.nnkpfile);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_method");
        item.annotation = "different implementation methods under Lcao basis set";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Integer";
        item.description = R"(Only available on LCAO basis, using different methods to generate "\.mmn" file and "\.amn" file.
* 1: Calculated using the lcao_in_pw method, the calculation accuracy can be improved by increasing ecutwfc to maintain consistency with the pw basis set results.
* 2: The overlap between atomic orbitals is calculated using grid integration. The radial grid points are generated using the Gauss-Legendre method, while the spherical grid points are generated using the Lebedev-Laikov method.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            /*
                       Developer's notes: on the repair of lcao_in_pw

                       lcao_in_pw is a special basis_type, for scf calculation,
                      it follows workflow of pw, but for nscf the toWannier90
                      calculation, the interface is in ESolver_KS_LCAO_elec,
                       therefore lcao_in_pw for towannier90 calculation follows
                      lcao.

                       In the future lcao_in_pw will have its own ESolver.

                       2023/12/22 use new psi_initializer to expand numerical
                      atomic orbitals, ykhuang
                   */
            if (para.input.towannier90 && para.input.basis_type == "lcao_in_pw")
            {
                para.input.wannier_method = 1;
            }
        };
        read_sync_int(input.wannier_method);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_spin");
        item.annotation = "calculate spin in wannier90 code interface";
        item.category = "Berry phase and wannier90 interface";
        item.type = "String";
        item.description = R"(The spin direction for the Wannier function calculation when nspin is set to 2
* up: Calculate spin up for the Wannier function.
* down: Calculate spin down for the Wannier function.)";
        item.default_value = "up";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.wannier_spin);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_mmn");
        item.annotation = "output .mmn file or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Write the "*.mmn" file or not.
* 0: don't write the "*.mmn" file.
* 1: write the "*.mmn" file.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wannier_mmn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_amn");
        item.annotation = "output .amn file or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Write the "*.amn" file or not.
* 0: don't write the "*.amn" file.
* 1: write the "*.amn" file.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wannier_amn);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_eig");
        item.annotation = "output .eig file or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Write the "*.eig" file or not.
* 0: don't write the "*.eig" file.
* 1: write the "*.eig" file.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wannier_eig);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_unk");
        item.annotation = "output UNK. file or not";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Write the "UNK.*" file or not.
* 0: don't write the "UNK.*" file.
* 1: write the "UNK.*" file.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wannier_unk);
        this->add_item(item);
    }
    {
        Input_Item item("out_wannier_wvfn_formatted");
        item.annotation = "output UNK. file in text format or in binary format";
        item.category = "Berry phase and wannier90 interface";
        item.type = "Boolean";
        item.description = R"(Write the "UNK.*" file in ASCII format or binary format.
* 0: write the "UNK.*" file in binary format.
* 1: write the "UNK.*" file in ASCII format (text file format).)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wannier_wvfn_formatted);
        this->add_item(item);
    }
}
} // namespace ModuleIO
