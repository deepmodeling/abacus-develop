#include "read_input.h"
#include "read_input_tool.h"


namespace ModuleIO
{
void ReadInput::item_sdft()
{
    {
        Input_Item item("method_sto");
        item.annotation = "1: slow and save memory, 2: fast and waste memory";
        read_sync_int(method_sto);
        this->add_item(item);
    }
    {
        Input_Item item("npart_sto");
        item.annotation = "Reduce memory when calculating Stochastic DOS";
        read_sync_int(npart_sto);
        this->add_item(item);
    }
    {
        Input_Item item("nbands_sto");
        item.annotation = "number of stochstic orbitals";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.nbndsto_str = strvalue;
            if(para.input.nbndsto_str != "all")
            {
                para.input.nbands_sto = std::stoi(para.input.nbndsto_str);
            }
        };
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            if(para.input.nbndsto_str == "all")
            {
                para.input.nbands_sto = 0;
            }
        };
        sync_string(nbndsto_str);
        add_int_bcast(nbands_sto); // Since "nbands_sto" has been assigned a value, it needs to be broadcasted
        this->add_item(item);
    }
    {
        Input_Item item("nche_sto");
        item.annotation = "Chebyshev expansion orders";
        read_sync_int(nche_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emin_sto");
        item.annotation = "trial energy to guess the lower bound of eigen energies of the Hamitonian operator";
        read_sync_double(emin_sto);
        this->add_item(item);
    }
    {
        Input_Item item("emax_sto");
        item.annotation = "trial energy to guess the upper bound of eigen energies of the Hamitonian operator";
        read_sync_double(emax_sto);
        this->add_item(item);
    }
    {
        Input_Item item("seed_sto");
        item.annotation = "the random seed to generate stochastic orbitals";
        read_sync_int(seed_sto);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_ecut");
        item.annotation = "maximum ecut to init stochastic bands";
        read_sync_double(initsto_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("initsto_freq");
        item.annotation = "frequency to generate new stochastic orbitals when running md";
        read_sync_int(initsto_freq);
        this->add_item(item);
    }
    {
        Input_Item item("cal_cond");
        item.annotation = "calculate electronic conductivities";
        read_sync_bool(cal_cond);
        this->add_item(item);
    }
    {
        Input_Item item("cond_che_thr");
        item.annotation = "control the error of Chebyshev expansions for conductivities";
        read_sync_double(cond_che_thr);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dw");
        item.annotation = "frequency interval for conductivities";
        read_sync_double(cond_dw);
        this->add_item(item);
    }
    {
        Input_Item item("cond_wcut");
        item.annotation = "cutoff frequency (omega) for conductivities";
        read_sync_double(cond_wcut);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dt");
        item.annotation = "t interval to integrate Onsager coefficiencies";
        read_sync_double(cond_dt);
        this->add_item(item);
    }
    {
        Input_Item item("cond_dtbatch");
        item.annotation = "exp(iH*dt*cond_dtbatch) is expanded with Chebyshev expansion";
        read_sync_int(cond_dtbatch);
        this->add_item(item);
    }
    {
        Input_Item item("cond_smear");
        item.annotation = "Smearing method for conductivities";
        read_sync_int(cond_smear);
        this->add_item(item);
    }
    {
        Input_Item item("cond_fwhm");
        item.annotation = "FWHM for conductivities";
        read_sync_double(cond_fwhm);
        this->add_item(item);
    }
    {
        Input_Item item("cond_nonlocal");
        item.annotation = "Nonlocal effects for conductivities";
        read_sync_bool(cond_nonlocal);
        this->add_item(item);
    }
}
} // namespace ModuleIO