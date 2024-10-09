#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_sdft()
{
    {
        Input_Item item("method_sto");
        item.annotation = "1: slow and save memory, 2: fast and waste memory";
        read_sync_int(input.method_sto);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.method_sto != 1 && para.input.method_sto != 2)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "method_sto should be 1 or 2");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("npart_sto");
        item.annotation = "Reduce memory when calculating Stochastic DOS";
        read_sync_int(input.npart_sto);
        this->add_item(item);
    }
    {
        Input_Item item("nbands_sto");
        item.annotation = "number of stochstic orbitals";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            std::string nbandsto_str = strvalue;
            if (nbandsto_str != "all")
            {
                para.input.nbands_sto = std::stoi(nbandsto_str);
            }
            else
            {
                para.input.nbands_sto = 0;
            }
        };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            // only do it when nbands_sto is set in INPUT
            if (item.is_read())
            {
                if (strvalue == "0" && para.input.esolver_type == "sdft")
                {
                    para.input.esolver_type = "ksdft";
                    ModuleBase::GlobalFunc::AUTO_SET("esolver_type", para.input.esolver_type);
                }
            }
        };
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.nbands_sto < 0 || para.input.nbands_sto > 100000)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "nbands_sto should be in the range of 0 to 100000");
            }
        };
        item.get_final_value = [](Input_Item& item, const Parameter& para) {
            if (item.str_values.size() == 0) // no nbands_sto in INPUT
            {
                item.final_value << para.input.nbands_sto;
            }
            else
            {
                item.final_value << item.str_values[0];
            }
        };
        add_int_bcast(input.nbands_sto);
        this->add_item(item);
    }
    {
        Input_Item item("nche_sto");
        item.annotation = "Chebyshev expansion orders";
        read_sync_int(input.nche_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emin_sto");
        item.annotation = "trial energy to guess the lower bound of eigen "
                          "energies of the Hamitonian operator";
        read_sync_double(input.emin_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emax_sto");
        item.annotation = "trial energy to guess the upper bound of eigen "
                          "energies of the Hamitonian operator";
        read_sync_double(input.emax_sto);
        this->add_item(item);
    }
    {
        Input_Item item("seed_sto");
        item.annotation = "the random seed to generate stochastic orbitals";
        read_sync_int(input.seed_sto);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_ecut");
        item.annotation = "maximum ecut to init stochastic bands";
        read_sync_double(input.initsto_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_freq");
        item.annotation = "frequency to generate new stochastic orbitals when running md";
        read_sync_int(input.initsto_freq);
        this->add_item(item);
    }
}
} // namespace ModuleIO