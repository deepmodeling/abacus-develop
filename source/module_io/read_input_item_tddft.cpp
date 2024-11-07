#include "module_base/constants.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_rt_tddft()
{ 
    // real time TDDFT
    {
        Input_Item item("td_force_dt");
        item.annotation = "time of force change";
        read_sync_double(input.td_force_dt);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext");
        item.annotation = "add extern potential or not";
        read_sync_bool(input.td_vext);
        this->add_item(item);
    }
    {
        Input_Item item("td_vext_dire");
        item.annotation = "extern potential direction";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_vext_dire = longstring(item.str_values);
        };
        sync_string(input.td_vext_dire);
        this->add_item(item);
    }
    {
        Input_Item item("init_vecpot_file");
        item.annotation = "init vector potential through file or not";
        read_sync_bool(input.init_vecpot_file);
        this->add_item(item);
    }
    {
        Input_Item item("td_print_eij");
        item.annotation = "print eij or not";
        read_sync_double(input.td_print_eij);
        this->add_item(item);
    }
    {
        Input_Item item("td_edm");
        item.annotation = "the method to calculate the energy density matrix";
        read_sync_int(input.td_edm);
        this->add_item(item);
    }
    {
        Input_Item item("td_propagator");
        item.annotation = "method of propagator";
        read_sync_int(input.propagator);
        this->add_item(item);
    }
    {
        Input_Item item("td_stype");
        item.annotation = "type of electric field in space domain";
        read_sync_int(input.td_stype);
        this->add_item(item);
    }
    {
        Input_Item item("td_ttype");
        item.annotation = "type of electric field in time domain";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_ttype = longstring(item.str_values);
        };
        sync_string(input.td_ttype);
        this->add_item(item);
    }
    {
        Input_Item item("td_tstart");
        item.annotation = " number of steps where electric field starts";
        read_sync_int(input.td_tstart);
        this->add_item(item);
    }
    {
        Input_Item item("td_tend");
        item.annotation = "number of steps where electric field ends";
        read_sync_int(input.td_tend);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut1");
        item.annotation = "cut1 of interval in length gauge";
        read_sync_double(input.td_lcut1);
        this->add_item(item);
    }
    {
        Input_Item item("td_lcut2");
        item.annotation = "cut2 of interval in length gauge";
        read_sync_double(input.td_lcut2);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_freq");
        item.annotation = "frequency (freq) of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_freq = longstring(item.str_values);
        };
        sync_string(input.td_gauss_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_phase");
        item.annotation = "phase of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_phase = longstring(item.str_values);
        };
        sync_string(input.td_gauss_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_sigma");
        item.annotation = "sigma of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_sigma = longstring(item.str_values);
        };
        sync_string(input.td_gauss_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_t0");
        item.annotation = "step number of time center (t0) of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_t0 = longstring(item.str_values);
        };
        sync_string(input.td_gauss_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_gauss_amp");
        item.annotation = "amplitude of Gauss type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_gauss_amp = longstring(item.str_values);
        };
        sync_string(input.td_gauss_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_freq");
        item.annotation = "frequency of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_freq = longstring(item.str_values);
        };
        sync_string(input.td_trape_freq);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_phase");
        item.annotation = "phase of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_phase = longstring(item.str_values);
        };
        sync_string(input.td_trape_phase);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t1");
        item.annotation = "t1 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t1 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t2");
        item.annotation = "t2 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t2 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_t3");
        item.annotation = "t3 of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_t3 = longstring(item.str_values);
        };
        sync_string(input.td_trape_t3);
        this->add_item(item);
    }
    {
        Input_Item item("td_trape_amp");
        item.annotation = "amplitude of Trapezoid type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trape_amp = longstring(item.str_values);
        };
        sync_string(input.td_trape_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq1");
        item.annotation = "frequency 1 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq1 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_freq1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_freq2");
        item.annotation = "frequency 2 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_freq2 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_freq2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase1");
        item.annotation = "phase 1 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase1 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_phase1);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_phase2");
        item.annotation = "phase 2 of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_phase2 = longstring(item.str_values);
        };
        sync_string(input.td_trigo_phase2);
        this->add_item(item);
    }
    {
        Input_Item item("td_trigo_amp");
        item.annotation = "amplitude of Trigonometric type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_trigo_amp = longstring(item.str_values);
        };
        sync_string(input.td_trigo_amp);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_t0");
        item.annotation = "t0 of Heaviside type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_t0 = longstring(item.str_values);
        };
        sync_string(input.td_heavi_t0);
        this->add_item(item);
    }
    {
        Input_Item item("td_heavi_amp");
        item.annotation = "amplitude of Heaviside type electric field";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.td_heavi_amp = longstring(item.str_values);
        };
        sync_string(input.td_heavi_amp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp");
        item.annotation = "change occupation or not";
        read_sync_bool(input.ocp);
        this->add_item(item);
    }
    {
        Input_Item item("ocp_set");
        item.annotation = "set occupation";
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
void ReadInput::item_lr_tddft()
{
    // Linear Responce TDDFT
    {
        Input_Item item("lr_nstates");
        item.annotation = "the number of 2-particle states to be solved";
        read_sync_int(input.lr_nstates);
        this->add_item(item);
    }
    {
        Input_Item item("nocc");
        item.annotation = "the number of occupied orbitals to form the 2-particle basis ( <= nelec/2)";
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
        read_sync_int(input.nvirt);
        this->add_item(item);
    }
    {
        Input_Item item("xc_kernel");
        item.annotation = "exchange correlation (XC) kernel for LR-TDDFT";
        read_sync_string(input.xc_kernel);
        this->add_item(item);
    }
    {
        Input_Item item("lr_init_xc_kernel");
        item.annotation = "The method to initalize the xc kernel";
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
        read_sync_string(input.lr_solver);
        this->add_item(item);
    }
    {
        Input_Item item("lr_thr");
        item.annotation = "convergence threshold of the LR-TDDFT eigensolver";
        read_sync_double(input.lr_thr);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lr");
        item.annotation = "whether to output the eigenvectors (excitation amplitudes) in the particle-hole basis";
        read_sync_bool(input.out_wfc_lr);
        this->add_item(item);
    }
    {
        Input_Item item("lr_unrestricted");
        item.annotation = "Whether to use unrestricted construction for LR-TDDFT";
        read_sync_bool(input.lr_unrestricted);
        this->add_item(item);
    }
    {
        Input_Item item("abs_wavelen_range");
        item.annotation = "the range of wavelength(nm) to output the absorption spectrum ";
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
        Input_Item item("abs_broadening");
        item.annotation = "the broadening (eta) for LR-TDDFT absorption spectrum";
        read_sync_double(input.abs_broadening);
        this->add_item(item);
    }
}
}