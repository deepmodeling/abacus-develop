#include "read_input.h"
#include "read_input_tool.h"


namespace ModuleIO
{
void ReadInput::item_postprocess()
{
    // 6. Smearing
    {
        Input_Item item("smearing_method");
        item.annotation = "type of smearing_method: gauss; fd; fixed; mp; mp2; mv";
        read_sync_string(smearing_method);
        this->add_item(item);
    }
    {
        Input_Item item("smearing_sigma");
        item.annotation = "energy range for smearing";
        read_sync_double(smearing_sigma);
        this->add_item(item);
    }

    // 7. Charge Mixing
    {
        Input_Item item("mixing_type");
        item.annotation = "plain; pulay; broyden";
        read_sync_string(mixing_mode);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta");
        item.annotation = "mixing parameter: 0 means no new charge";
        read_sync_double(mixing_beta);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_ndim");
        item.annotation = "mixing dimension in pulay or broyden";
        read_sync_int(mixing_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_restart");
        item.annotation = "threshold to restart mixing during SCF";
        read_sync_double(mixing_restart);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0");
        item.annotation = "mixing parameter in kerker";
        read_sync_double(mixing_gg0);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta_mag");
        item.annotation = "mixing parameter for magnetic density";
        read_sync_double(mixing_beta_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_mag");
        item.annotation = "mixing parameter in kerker";
        read_sync_double(mixing_gg0_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_min");
        item.annotation = "the minimum kerker coefficient";
        read_sync_double(mixing_gg0_min);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_angle");
        item.annotation = "angle mixing parameter for non-colinear calculations";
        read_sync_double(mixing_angle);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_tau");
        item.annotation = "whether to mix tau in mGGA calculation";
        read_sync_bool(mixing_tau);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dftu");
        item.annotation = "whether to mix locale in DFT+U calculation";
        read_sync_bool(mixing_dftu);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dmr");
        item.annotation = "whether to mix real-space density matrix";
        read_sync_bool(mixing_dmr);
        this->add_item(item);
    }

    // 8. DOS
    {
        Input_Item item("dos_emin_ev");
        item.annotation = "minimal range for dos";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emin_ev = doublevalue;
            para.input.dos_setemin = true;
        };
        sync_double(dos_emin_ev);
        add_bool_bcast(dos_setemin); //Since "dos_setemin" has been assigned a value, it needs to be broadcasted
        this->add_item(item);
    }
    {
        Input_Item item("dos_emax_ev");
        item.annotation = "maximal range for dos";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emax_ev = doublevalue;
            para.input.dos_setemax = true;
        };
        sync_double(dos_emax_ev);
        add_bool_bcast(dos_setemax); //Since "dos_setemax" has been assigned a value, it needs to be broadcasted
        this->add_item(item);
    }
    {
        Input_Item item("dos_edelta_ev");
        item.annotation = "delta energy for dos";
        read_sync_double(dos_edelta_ev);
        this->add_item(item);
    }
    {
        Input_Item item("dos_scale");
        item.annotation = "scale dos range by";
        read_sync_double(dos_scale);
        this->add_item(item);
    }
    {
        Input_Item item("dos_sigma");
        item.annotation = "gauss b coefficeinet(default=0.07)";
        read_sync_double(dos_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("dos_nche");
        item.annotation = "orders of Chebyshev expansions for dos";
        read_sync_int(dos_nche);
        this->add_item(item);
    }
}
} // namespace ModuleIO