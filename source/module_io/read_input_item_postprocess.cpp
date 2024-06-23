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
        read_bcast_string(smearing_method);
        this->add_item(item);
    }
    {
        Input_Item item("smearing_sigma");
        item.annotation = "energy range for smearing";
        read_bcast_double(smearing_sigma);
        this->add_item(item);
    }

    // 7. Charge Mixing
    {
        Input_Item item("mixing_type");
        item.annotation = "plain; pulay; broyden";
        read_bcast_string(mixing_mode);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta");
        item.annotation = "mixing parameter: 0 means no new charge";
        read_bcast_double(mixing_beta);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_ndim");
        item.annotation = "mixing dimension in pulay or broyden";
        read_bcast_int(mixing_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_restart");
        item.annotation = "threshold to restart mixing during SCF";
        read_bcast_double(mixing_restart);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0");
        item.annotation = "mixing parameter in kerker";
        read_bcast_double(mixing_gg0);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_beta_mag");
        item.annotation = "mixing parameter for magnetic density";
        read_bcast_double(mixing_beta_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_mag");
        item.annotation = "mixing parameter in kerker";
        read_bcast_double(mixing_gg0_mag);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_gg0_min");
        item.annotation = "the minimum kerker coefficient";
        read_bcast_double(mixing_gg0_min);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_angle");
        item.annotation = "angle mixing parameter for non-colinear calculations";
        read_bcast_double(mixing_angle);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_tau");
        item.annotation = "whether to mix tau in mGGA calculation";
        read_bcast_bool(mixing_tau);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dftu");
        item.annotation = "whether to mix locale in DFT+U calculation";
        read_bcast_bool(mixing_dftu);
        this->add_item(item);
    }
    {
        Input_Item item("mixing_dmr");
        item.annotation = "whether to mix real-space density matrix";
        read_bcast_bool(mixing_dmr);
        this->add_item(item);
    }

    // 8. DOS
    {
        Input_Item item("dos_emin_ev");
        item.annotation = "minimal range for dos";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.dos_emin_ev = doublevalue;
            para.dos_setemin = true;
        };
        bcast_double_param(dos_emin_ev);
        bcast_bool_param(dos_setemin);
        this->add_item(item);
    }
    {
        Input_Item item("dos_emax_ev");
        item.annotation = "maximal range for dos";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.dos_emax_ev = doublevalue;
            para.dos_setemax = true;
        };
        bcast_double_param(dos_emax_ev);
        bcast_bool_param(dos_setemax);
        this->add_item(item);
    }
    {
        Input_Item item("dos_edelta_ev");
        item.annotation = "delta energy for dos";
        read_bcast_double(dos_edelta_ev);
        this->add_item(item);
    }
    {
        Input_Item item("dos_scale");
        item.annotation = "scale dos range by";
        read_bcast_double(dos_scale);
        this->add_item(item);
    }
    {
        Input_Item item("dos_sigma");
        item.annotation = "gauss b coefficeinet(default=0.07)";
        read_bcast_double(dos_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("dos_nche");
        item.annotation = "orders of Chebyshev expansions for dos";
        read_bcast_int(dos_nche);
        this->add_item(item);
    }
}
} // namespace ModuleIO