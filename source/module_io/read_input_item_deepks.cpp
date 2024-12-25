#include "module_base/constants.h"
#include "module_parameter/parameter.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"
namespace ModuleIO
{
void ReadInput::item_deepks()
{
    {
        Input_Item item("deepks_out_labels");
        item.annotation = ">0 compute descriptor for deepks";
        read_sync_bool(input.deepks_out_labels);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_scf");
        item.annotation = ">0 add V_delta to Hamiltonian";
        read_sync_bool(input.deepks_scf);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
        #ifndef __DEEPKS
            if (PARAM.inp.deepks_scf || PARAM.inp.deepks_out_labels || 
                PARAM.inp.deepks_bandgap || PARAM.inp.deepks_v_delta)
            {
                ModuleBase::WARNING_QUIT("Input_conv", "please compile with DeePKS");
            }
        #endif
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_equiv");
        item.annotation = "whether to use equivariant version of DeePKS";
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
        Input_Item item("deepks_bandgap");
        item.annotation = ">0 for bandgap label";
        read_sync_bool(input.deepks_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_v_delta");
        item.annotation = ">0 for v_delta label. when output, 1 for v_delta_precalc, 2 for phialpha and grad_evdm ( can save memory )";
        read_sync_int(input.deepks_v_delta);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_unittest");
        item.annotation = "if set 1, prints intermediate quantities that shall "
                          "be used for making unit test";
        read_sync_bool(input.deepks_out_unittest);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.deepks_out_unittest)
            {
                para.input.deepks_out_labels = true;
                para.input.deepks_scf = true;

            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.deepks_out_unittest){
                if (para.input.cal_force != 1 || para.input.cal_stress != 1)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "force and stress are required in generating deepks unittest");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("deepks_model");
        item.annotation = "file dir of traced pytorch model: 'model.ptg";
        read_sync_string(input.deepks_model);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_lmax");
        item.annotation = "lmax used in generating spherical bessel functions";
        read_sync_int(input.bessel_descriptor_lmax);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_ecut");
        item.annotation = "energy cutoff for spherical bessel functions(Ry)";
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
        read_sync_double(input.bessel_descriptor_tolerence);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_rcut");
        item.annotation = "radial cutoff for spherical bessel functions(a.u.)";
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
        read_sync_bool(input.bessel_descriptor_smooth);
        this->add_item(item);
    }
    {
        Input_Item item("bessel_descriptor_sigma");
        item.annotation = "sphereical bessel smearing_sigma";
        read_sync_double(input.bessel_descriptor_sigma);
        this->add_item(item);
    }
}
} // namespace ModuleIO