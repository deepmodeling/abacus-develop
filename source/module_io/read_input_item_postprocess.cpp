#include "module_base/tool_quit.h"
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
    {
        // Energy range for smearing,
        //`smearing_sigma` = 1/2 *kB* `smearing_sigma_temp`.
        Input_Item tmp_item("smearing_sigma_temp");
        tmp_item.read_value
            = [](const Input_Item& item, Parameter& para) { para.input.smearing_sigma = 3.166815e-6 * doublevalue; };
        // only to set smearing_sigma, so no need to write to output INPUT file
        // or bcast.
        this->add_item(tmp_item);
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
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mixing_beta < 0.0)
            {
                if (para.input.nspin == 1)
                {
                    para.input.mixing_beta = 0.8;
                }
                else if (para.input.nspin == 2)
                {
                    para.input.mixing_beta = 0.4;
                    para.input.mixing_beta_mag = 1.6;
                    para.input.mixing_gg0_mag = 0.0;
                }
                else if (para.input.nspin == 4) // I will add this
                {
                    para.input.mixing_beta = 0.4;
                    para.input.mixing_beta_mag = 1.6;
                    para.input.mixing_gg0_mag = 0.0;
                }
            }
        };
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
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.mixing_beta_mag < 0.0)
            {
                if (para.input.nspin == 2 || para.input.nspin == 4)
                {
                    if (para.input.mixing_beta <= 0.4)
                    {
                        para.input.mixing_beta_mag = 4 * para.input.mixing_beta;
                    }
                    else
                    {
                        para.input.mixing_beta_mag = 1.6; // 1.6 can be discussed
                    }
                }
            }
        };
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
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emin_ev = doublevalue;
            para.input.sup.dos_setemin = true;
        };
        sync_double(dos_emin_ev);
        add_bool_bcast(sup.dos_setemin);
        this->add_item(item);
    }
    {
        Input_Item item("dos_emax_ev");
        item.annotation = "maximal range for dos";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            para.input.dos_emax_ev = doublevalue;
            para.input.sup.dos_setemax = true;
        };
        sync_double(dos_emax_ev);
        add_bool_bcast(sup.dos_setemax);
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