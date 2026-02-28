#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_sdft()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("method_sto");
        item.annotation = "1: slow and save memory, 2: fast and waste memory";
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer";
        item.description = R"(Different methods to do stochastic DFT
* 1: Calculate twice, this method cost less memory but is slower.
* 2: Calculate once but needs much more memory. This method is much faster. Besides, it calculates with a smaller nche_sto. However, when the memory is not enough, only method 1 can be used.
* other: use 2)";
        item.default_value = "2";
        item.unit = "";
        item.availability = "esolver_type = sdft";
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
        Input_Item item("nbands_sto");
        item.annotation = "number of stochstic orbitals";
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer or string";
        item.description = R"(The number of stochastic orbitals
* > 0: Perform stochastic DFT. Increasing the number of bands improves accuracy and reduces stochastic errors; To perform mixed stochastic-deterministic DFT, you should set nbands, which represents the number of KS orbitals.
* 0: Perform Kohn-Sham DFT.
* all: All complete basis sets are used to replace stochastic orbitals with the Chebyshev method (CT), resulting in the same results as KSDFT without stochastic errors.)";
        item.default_value = "256";
        item.unit = "";
        item.availability = "esolver_type = sdft";
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
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer";
        item.description = "Chebyshev expansion orders for stochastic DFT.";
        item.default_value = "100";
        item.unit = "";
        item.availability = "esolver_type = sdft";
        read_sync_int(input.nche_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emin_sto");
        item.annotation = "trial energy to guess the lower bound of eigen "
                          "energies of the Hamitonian operator";
        item.category = "Electronic structure (SDFT)";
        item.type = "Real";
        item.description = "Trial energy to guess the lower bound of eigen energies of the Hamiltonian Operator.";
        item.default_value = "0.0";
        item.unit = "Ry";
        item.availability = "esolver_type = sdft";
        read_sync_double(input.emin_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emax_sto");
        item.annotation = "trial energy to guess the upper bound of eigen "
                          "energies of the Hamitonian operator";
        item.category = "Electronic structure (SDFT)";
        item.type = "Real";
        item.description = "Trial energy to guess the upper bound of eigen energies of the Hamiltonian Operator.";
        item.default_value = "0.0";
        item.unit = "Ry";
        item.availability = "esolver_type = sdft";
        read_sync_double(input.emax_sto);
        this->add_item(item);
    }
    {
        Input_Item item("seed_sto");
        item.annotation = "the random seed to generate stochastic orbitals";
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer";
        item.description = R"(The random seed to generate stochastic orbitals.
* >= 0: Stochastic orbitals have the form of exp(i*theta), where theta is a uniform distribution in [0, 2*pi).
* 0: the seed is decided by time(NULL).
* <= -1: Stochastic orbitals have the form of +1 or -1 with equal probability.
* -1: the seed is decided by time(NULL).)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "esolver_type = sdft";
        read_sync_int(input.seed_sto);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_ecut");
        item.annotation = "maximum ecut to init stochastic bands";
        item.category = "Electronic structure (SDFT)";
        item.type = "Real";
        item.description = R"(Stochastic wave functions are initialized in a large box generated by "4*initsto_ecut". initsto_ecut should be larger than ecutwfc. In this method, SDFT results are the same when using different cores. Besides, coefficients of the same G are the same when ecutwfc is rising to initsto_ecut. If it is smaller than ecutwfc, it will be turned off.)";
        item.default_value = "0.0";
        item.unit = "Ry";
        item.availability = "esolver_type = sdft";
        read_sync_double(input.initsto_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_freq");
        item.annotation = "frequency to generate new stochastic orbitals when running md";
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer";
        item.description = R"(Frequency (once each initsto_freq steps) to generate new stochastic orbitals when running md.
* positive integer: Update stochastic orbitals
* 0: Never change stochastic orbitals.)";
        item.default_value = "0";
        item.unit = "";
        item.availability = "esolver_type = sdft";
        read_sync_int(input.initsto_freq);
        this->add_item(item);
    }
    {
        Input_Item item("npart_sto");
        item.annotation = "Reduce memory when calculating Stochastic DOS";
        item.category = "Electronic structure (SDFT)";
        item.type = "Integer";
        item.description = "Make memory cost to 1/npart_sto times of the previous one when running the post process of SDFT like DOS or conductivities.";
        item.default_value = "1";
        item.unit = "";
        item.availability = "method_sto = 2 and out_dos = 1 or cal_cond = True";
        read_sync_int(input.npart_sto);
        this->add_item(item);
    }
}
} // namespace ModuleIO
