#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_rt_tddft()
{ 
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("estep_per_md");
        item.annotation = "steps of force change";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = "The number of electronic propagation steps between two ionic steps.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.estep_per_md);
        this->add_item(item);
    }
    // real time TDDFT
    {
        Input_Item item("td_dt");
        item.annotation = "time step for evolving wavefunction";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = "The time step used in electronic propagation. Setting td_dt will reset the value of md_dt to td_dt * estep_per_md.";
        item.default_value = "md_dt / estep_per_md";
        item.unit = "fs";
        item.availability = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.td_dt == -1.0)
            {
                GlobalV::ofs_running << "td_dt don't exist, set td_dt with md_dt" << std::endl;
                para.input.td_dt = para.input.mdp.md_dt / para.input.estep_per_md;
            }
        };
        read_sync_double(input.td_dt);
        this->add_item(item);
    }
    {
        Input_Item item("td_edm");
        item.annotation = "the method to calculate the energy density matrix";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Method to calculate the energy-density matrix, mainly affects the calculation of force and stress.
* 0: Using the original formula.
* 1: Using the formula for ground state (deprecated). Note that this usually does not hold if wave function is not the eigenstate of the Hamiltonian.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.td_edm);
        this->add_item(item);
    }
    {
        Input_Item item("td_print_eij");
        item.annotation = "print eij or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = R"(Controls the printing of Hamiltonian matrix elements.
* < 0: Suppress all output.
* >= 0: Print only elements with either i or j exceeding td_print_eij.)";
        item.default_value = "-1";
        item.unit = "Ry";
        item.availability = "";
        read_sync_double(input.td_print_eij);
        this->add_item(item);
    }
    {
        Input_Item item("td_propagator");
        item.annotation = "method of propagator";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Methods of electronic propagation.
* 0: Crank-Nicolson, based on matrix inversion.
* 1: 4th-order Taylor expansion of exponential.
* 2: Enforced time-reversal symmetry (ETRS).
* 3: Crank-Nicolson, based on solving linear equation.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.propagator);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext");
        item.annotation = "add extern potential or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(* True: Add a laser-material interaction (external electric field).
* False: No external electric field.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.td_vext);
        this->add_item(item);
    }
    // {
    //     Input_Item item("td_vext_dire");
    //     item.annotation = "extern potential direction";
    //     item.read_value = [](const Input_Item& item, Parameter& para) {
    //         para.input.td_vext_dire = longstring(item.str_values);
    //     };
    //     sync_string(input.td_vext_dire);
    //     this->add_item(item);
    // }
    {
        Input_Item item("td_vext_dire");
        item.annotation = "extern potential direction";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = R"(Specifies the direction(s) of the external electric field when td_vext is enabled. For example, td_vext_dire 1 2 indicates that external electric fields are applied to both the x and y directions simultaneously. Electric field parameters can also be written as strings. For example, td_gauss_phase 0 1.5707963 indicates that the Gaussian type electric fields in the x and y directions have a phase delay of pi/2.
* 1: The external field direction is along the x-axis.
* 2: The external field direction is along the y-axis.
* 3: The external field direction is along the z-axis.)";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.td_vext_dire);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_intvec_bcast(input.td_vext_dire, para.input.td_vext_dire.size(), 0);
        this->add_item(item);
    }
    {
        Input_Item item("td_stype");
        item.annotation = "type of electric field in space domain";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = R"(Type of electric field in the space domain, i.e. the gauge of the electric field.
* 0: Length gauge.
* 1: Velocity gauge.
* 2: Hybrid gauge. See J. Chem. Theory Comput. 2025, 21, 3335-3341 for more information.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.td_stype);
        this->add_item(item);
    }
    {
        Input_Item item("td_ttype");
        item.annotation = "type of electric field in time domain";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = R"(Type of electric field in the time domain.
* 0: Gaussian type function.
* 1: Trapezoid type function.
* 2: Trigonometric type function.
* 3: Heaviside type function.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_ttype = longstring(item.str_values);
        };
        sync_string(input.td_ttype);
        this->add_item(item);
    }
    {
        Input_Item item("td_tstart");
        item.annotation = " number of steps where electric field starts";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = "The initial time step when the time-dependent electric field is activated.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.td_tstart);
        this->add_item(item);
    }
    {
        Input_Item item("td_tend");
        item.annotation = "number of steps where electric field ends";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Integer";
        item.description = "The final time step when the time-dependent electric field is deactivated. The field remains active between td_tstart and td_tend.";
        item.default_value = "1000";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.td_tend);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut1");
        item.annotation = "cut1 of interval in length gauge";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = "The lower bound of the interval in the length gauge RT-TDDFT, where the coordinate is the fractional coordinate.";
        item.default_value = "0.05";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.td_lcut1);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut2");
        item.annotation = "cut2 of interval in length gauge";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Real";
        item.description = "The upper bound of the interval in the length gauge RT-TDDFT, where the coordinate is the fractional coordinate.";
        item.default_value = "0.95";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.td_lcut2);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_freq");
        item.annotation = "frequency (freq) of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Frequency of the Gaussian type electric field.";
        item.default_value = "22.13";
        item.unit = "1/fs";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_freq = longstring(item.str_values);
        };
        sync_string(input.td_gauss_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_phase");
        item.annotation = "phase of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Phase of the Gaussian type electric field.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_phase = longstring(item.str_values);
        };
        sync_string(input.td_gauss_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_sigma");
        item.annotation = "sigma of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Pulse width (standard deviation) of the Gaussian type electric field.";
        item.default_value = "30.0";
        item.unit = "fs";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_sigma = longstring(item.str_values);
        };
        sync_string(input.td_gauss_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_t0");
        item.annotation = "step number of time center (t0) of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Step number of the time center of the Gaussian type electric field.";
        item.default_value = "100";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_t0 = longstring(item.str_values);
        };
        sync_string(input.td_gauss_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_amp");
        item.annotation = "amplitude of Gauss type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Amplitude of the Gaussian type electric field.";
        item.default_value = "0.25";
        item.unit = "V/Angstrom";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_amp = longstring(item.str_values);
        };
        sync_string(input.td_gauss_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_freq");
        item.annotation = "frequency of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Frequency of the trapezoid type electric field.";
        item.default_value = "1.60";
        item.unit = "1/fs";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_freq = longstring(item.str_values);
        };
        sync_string(input.td_trape_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_phase");
        item.annotation = "phase of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Phase of the trapezoid type electric field.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_phase = longstring(item.str_values);
        };
        sync_string(input.td_trape_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t1");
        item.annotation = "t1 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Step number of the time interval t1 of the trapezoid type electric field.";
        item.default_value = "1875";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t1 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t2");
        item.annotation = "t2 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Step number of the time interval t2 of the trapezoid type electric field.";
        item.default_value = "5625";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t2 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t3");
        item.annotation = "t3 of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Step number of the time interval t3 of the trapezoid type electric field.";
        item.default_value = "7500";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t3 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t3);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_amp");
        item.annotation = "amplitude of Trapezoid type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Amplitude of the trapezoid type electric field.";
        item.default_value = "2.74";
        item.unit = "V/Angstrom";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_amp = longstring(item.str_values);
        };
        sync_string(input.td_trape_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq1");
        item.annotation = "frequency 1 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Frequency 1 of the trigonometric type electric field.";
        item.default_value = "1.164656";
        item.unit = "1/fs";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq1 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_freq1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq2");
        item.annotation = "frequency 2 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Frequency 2 of the trigonometric type electric field.";
        item.default_value = "0.029116";
        item.unit = "1/fs";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq2 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_freq2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase1");
        item.annotation = "phase 1 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Phase 1 of the trigonometric type electric field.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase1 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_phase1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase2");
        item.annotation = "phase 2 of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Phase 2 of the trigonometric type electric field.";
        item.default_value = "0.0";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase2 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_phase2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_amp");
        item.annotation = "amplitude of Trigonometric type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Amplitude of the trigonometric type electric field.";
        item.default_value = "2.74";
        item.unit = "V/Angstrom";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_amp = longstring(item.str_values);
        };
        sync_string(input.td_trigo_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_t0");
        item.annotation = "t0 of Heaviside type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Step number of the switch time of the Heaviside type electric field.";
        item.default_value = "100";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_t0 = longstring(item.str_values);
        };
        sync_string(input.td_heavi_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_amp");
        item.annotation = "amplitude of Heaviside type electric field";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = "Amplitude of the Heaviside type electric field.";
        item.default_value = "1.0";
        item.unit = "V/Angstrom";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_amp = longstring(item.str_values);
        };
        sync_string(input.td_heavi_amp);
        this->add_item(item);
    }
    {
        Input_Item item("init_vecpot_file");
        item.annotation = "init vector potential through file or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(Initialize vector potential through file or not.
* True: Initialize vector potential from file At.dat (unit: a.u.). It consists of four columns, representing the step number and vector potential on each direction.
* False: Calculate vector potential by integrating the electric field.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.init_vecpot_file);
        this->add_item(item);
    }
    {
        Input_Item item("ocp");
        item.annotation = "change occupation or not";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "Boolean";
        item.description = R"(* True: Fixes the band occupations based on the values specified in ocp_set.
* False: Does not fix the band occupations.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.ocp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp_set");
        item.annotation = "set occupation";
        item.category = "RT-TDDFT: Real-Time Time-Dependent Density Functional Theory";
        item.type = "String";
        item.description = R"(If ocp is set to 1, ocp_set must be provided as a string specifying the occupation numbers for each band across all k-points. The format follows a space-separated pattern, where occupations are assigned sequentially to bands for each k-point. A shorthand notation Nx can be used to repeat a value x for N bands.
* Example:
  1 10*1 0 1 represents occupations for 13 bands, where the 12th band is fully unoccupied (0), and all others are occupied (1).
* For a system with multiple k-points, the occupations must be specified for all k-points, following their order in the output file kpoints (may lead to fractional occupations).
* Incorrect specification of ocp_set could lead to inconsistencies in electron counting, causing the calculation to terminate with an error.)";
        item.default_value = "None";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.ocp_kb);
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if(item.is_read())
            {
                item.final_value.str(longstring(item.str_values));
            }
        };
        add_doublevec_bcast(input.ocp_kb, para.input.ocp_kb.size(), 0.0);
        this->add_item(item);
    }


}
void ReadInput::item_tdofdft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    // TD-OFDFT
    {
        Input_Item item("of_cd");
        item.annotation = "add CD Potential or not";
        item.category = "TDOFDFT: time dependent orbital free density functional theory";
        item.type = "Boolean";
        item.description = R"(Added the current dependent(CD) potential. (https://doi.org/10.1103/PhysRevB.98.144302)
* True: Added the CD potential.
* False: Not added the CD potential.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "TDOFDFT";
        read_sync_bool(input.of_cd);
        this->add_item(item);
    }
    {
        Input_Item item("of_mcd_alpha");
        item.annotation = "parameter of modified CD Potential";
        item.category = "TDOFDFT: time dependent orbital free density functional theory";
        item.type = "Real";
        item.description = "The value of the parameter alpha in modified CD potential method. mCDPotential=alpha*CDPotential (proposed in paper PhysRevB.98.144302)";
        item.default_value = "1.0";
        item.unit = "";
        item.availability = "TDOFDFT";
        read_sync_double(input.of_mCD_alpha);
        this->add_item(item);
    }
}
void ReadInput::item_lr_tddft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("xc_kernel");
        item.annotation = "exchange correlation (XC) kernel for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = "The exchange-correlation kernel used in the calculation. Currently supported: RPA, LDA, PBE, HSE, HF.";
        item.default_value = "LDA";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.xc_kernel);
        this->add_item(item);
    }
    {
        Input_Item item("lr_init_xc_kernel");
        item.annotation = "The method to initalize the xc kernel";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Vector of String (>=1 values)";
        item.description = R"(The method to initalize the xc kernel.
* "default": Calculate xc kernel from the ground-state charge density.
* "file": Read the xc kernel on grid from the provided files. The following words should be the paths of ".cube" files, where the first 1 (nspin==1) or 3 (nspin==2, namely spin-aa, spin-ab and spin-bb) will be read in. The parameter xc_kernel will be invalid. Now only LDA-type kernel is supported as the potential will be calculated by directly multiplying the transition density.
* "from_charge_file": Calculate fxc from the charge density read from the provided files. The following words should be the paths of ".cube" files, where the first nspin files will be read in.)";
        item.default_value = "\"default\"";
        item.unit = "";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            auto& ifxc = para.input.lr_init_xc_kernel;
            for (int i = 0; i < count; i++) { ifxc.push_back(item.str_values[i]); }
            };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.lr_init_xc_kernel.empty()) { para.input.lr_init_xc_kernel.push_back("default"); }
            };
        sync_stringvec(input.lr_init_xc_kernel, para.input.lr_init_xc_kernel.size(), "default");
        this->add_item(item);
    }
    {
        Input_Item item("lr_solver");
        item.annotation = "the eigensolver for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = R"(The method to solve the Casida equation in LR-TDDFT under Tamm-Dancoff approximation (TDA).
* dav/dav_subspace/cg: Construct and diagonalize the Hamiltonian matrix iteratively with Davidson/Non-ortho-Davidson/CG algorithm.
* lapack: Construct the full matrix and directly diagonalize with LAPACK.
* spectrum: Calculate absorption spectrum only without solving Casida equation.)";
        item.default_value = "dav";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.lr_solver);
        this->add_item(item);
    }
    {
        Input_Item item("lr_thr");
        item.annotation = "convergence threshold of the LR-TDDFT eigensolver";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real";
        item.description = "The convergence threshold of iterative diagonalization solver for LR-TDDFT. It is a pure-math number with the same meaning as pw_diag_thr, but since the Casida equation is a one-shot eigenvalue problem, it is also the convergence threshold of LR-TDDFT.";
        item.default_value = "1e-2";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.lr_thr);
        this->add_item(item);
    }
    {
        Input_Item item("nocc");
        item.annotation = "the number of occupied orbitals to form the 2-particle basis ( <= nelec/2)";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = R"(The number of occupied orbitals (up to HOMO) used in the LR-TDDFT calculation.
* Note: If the value is illegal ( > nelec/2 or <= 0), it will be autoset to nelec/2.)";
        item.default_value = "nband";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.nocc);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const int nocc_default = std::max(static_cast<int>(para.input.nelec + 1) / 2, para.input.nbands);
            if (para.input.nocc <= 0 || para.input.nocc > nocc_default) { para.input.nocc = nocc_default; }
            };
        this->add_item(item);
    }
    {
        Input_Item item("nvirt");
        item.annotation = "the number of virtual orbitals to form the 2-particle basis (nocc + nvirt <= nbands)";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = "The number of virtual orbitals (starting from LUMO) used in the LR-TDDFT calculation.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.nvirt);
        this->add_item(item);
    }
    // Linear Responce TDDFT
    {
        Input_Item item("lr_nstates");
        item.annotation = "the number of 2-particle states to be solved";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Integer";
        item.description = "The number of 2-particle states to be solved.";
        item.default_value = "0";
        item.unit = "";
        item.availability = "";
        read_sync_int(input.lr_nstates);
        this->add_item(item);
    }
    {
        Input_Item item("lr_unrestricted");
        item.annotation = "Whether to use unrestricted construction for LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Boolean";
        item.description = R"(Whether to use unrestricted construction for LR-TDDFT (the matrix size will be doubled).
* True: Always use unrestricted LR-TDDFT.
* False: Use unrestricted LR-TDDFT only when the system is open-shell.)";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.lr_unrestricted);
        this->add_item(item);
    }
    {
        Input_Item item("abs_wavelen_range");
        item.annotation = "the range of wavelength(nm) to output the absorption spectrum ";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real Real";
        item.description = "The range of the wavelength for the absorption spectrum calculation.";
        item.default_value = "0.0 0.0";
        item.unit = "nm";
        item.availability = "";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            for (int i = 0; i < count; i++)
            {
                para.input.abs_wavelen_range.push_back(std::stod(item.str_values[i]));
            }
            };
        sync_doublevec(input.abs_wavelen_range, 2, 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lr");
        item.annotation = "whether to output the eigenvectors (excitation amplitudes) in the particle-hole basis";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Boolean";
        item.description = "Whether to output the eigenstates (excitation energy) and eigenvectors (excitation amplitude) of the LR-TDDFT calculation. The output files are OUT.{suffix}/Excitation_Amplitude_${processor_rank}.dat.";
        item.default_value = "False";
        item.unit = "";
        item.availability = "";
        read_sync_bool(input.out_wfc_lr);
        this->add_item(item);
    }
    {
        Input_Item item("abs_gauge");
        item.annotation = "whether to use length or velocity gauge to calculate the absorption spectrum in LR-TDDFT";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "String";
        item.description = "Whether to use length or velocity gauge to calculate the absorption spectrum in LR-TDDFT.";
        item.default_value = "length";
        item.unit = "";
        item.availability = "";
        read_sync_string(input.abs_gauge);
        this->add_item(item);
    }
    {
        Input_Item item("abs_broadening");
        item.annotation = "the broadening (eta) for LR-TDDFT absorption spectrum";
        item.category = "Linear Response TDDFT (Under Development Feature)";
        item.type = "Real";
        item.description = "The broadening factor for the absorption spectrum calculation.";
        item.default_value = "0.01";
        item.unit = "";
        item.availability = "";
        read_sync_double(input.abs_broadening);
        this->add_item(item);
    }
}
}
