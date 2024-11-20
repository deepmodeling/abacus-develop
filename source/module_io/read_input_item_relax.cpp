#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{
void ReadInput::item_relax()
{
    {
        Input_Item item("relax_method");
        item.annotation = "cg; bfgs; sd; cg; cg_bfgs;";
        read_sync_string(input.relax_method);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            const std::vector<std::string> relax_methods = {"cg", "bfgs", "sd", "cg_bfgs"};
            if (std::find(relax_methods.begin(), relax_methods.end(), para.input.relax_method) == relax_methods.end())
            {
                const std::string warningstr = nofound_str(relax_methods, "relax_method");
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("relax_new");
        item.annotation = "whether to use the new relaxation method";
        read_sync_bool(input.relax_new);
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.relax_new && para.input.relax_method != "cg")
            {
                para.input.relax_new = false;
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("relax");
        item.annotation = "allow relaxation along the specific direction";
        read_sync_bool(input.relax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_scale_force");
        item.annotation = "controls the size of the first CG step if relax_new is true";
        read_sync_double(input.relax_scale_force);
        this->add_item(item);
    }
    {
        Input_Item item("relax_nmax");
        item.annotation = "number of ion iteration steps";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::string& calculation = para.input.calculation;
            const std::vector<std::string> singlelist
                = {"scf", "nscf", "get_S", "get_pchg", "get_wf", "test_memory", "test_neighbour", "gen_bessel"};
            if (std::find(singlelist.begin(), singlelist.end(), calculation) != singlelist.end())
            {
                para.input.relax_nmax = 1;
            }
            else if (calculation == "relax" || calculation == "cell-relax")
            {
                if (para.input.relax_nmax == 0) // default value
                {
                    para.input.relax_nmax = 50;
                }
            }
        };
        read_sync_int(input.relax_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_cg_thr");
        item.annotation = "threshold for switching from cg to bfgs, unit: eV/Angstrom";
        read_sync_double(input.relax_cg_thr);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr");
        item.annotation = "force threshold, unit: Ry/Bohr";
        // read_sync_double(input.force_thr);
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.force_thr = doublevalue; };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.force_thr == -1 && para.input.force_thr_ev == -1)
            {
                para.input.force_thr = 1.0e-3; // default value
                para.input.force_thr_ev = para.input.force_thr * 13.6058 / 0.529177;
            }
            else if (para.input.force_thr == -1 && para.input.force_thr_ev != -1)
            {
                para.input.force_thr = para.input.force_thr_ev / 13.6058 * 0.529177;
            }
            else
            {
                // if both force_thr and force_thr_ev are set, use force_thr
                ModuleBase::WARNING("ReadInput", "both force_thr and force_thr_ev are set, use force_thr");
                para.input.force_thr_ev = para.input.force_thr * 13.6058 / 0.529177;
            }
        };
        sync_double(input.force_thr);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev");
        item.annotation = "force threshold, unit: eV/Angstrom";
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.force_thr_ev = doublevalue; };
        sync_double(input.force_thr_ev);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev2");
        item.annotation = "force invalid threshold, unit: eV/Angstrom";
        read_sync_double(input.force_thr_ev2);
        this->add_item(item);
    }
    {
        Input_Item item("stress_thr");
        item.annotation = "stress threshold";
        read_sync_double(input.stress_thr);
        this->add_item(item);
    }
    {
        Input_Item item("press1");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press1);
        this->add_item(item);
    }
    {
        Input_Item item("press2");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press2);
        this->add_item(item);
    }
    {
        Input_Item item("press3");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(input.press3);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w1");
        item.annotation = "wolfe condition 1 for bfgs";
        read_sync_double(input.relax_bfgs_w1);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w2");
        item.annotation = "wolfe condition 2 for bfgs";
        read_sync_double(input.relax_bfgs_w2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmax");
        item.annotation = "maximal trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmin");
        item.annotation = "minimal trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_rmin);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_init");
        item.annotation = "initial trust radius, unit: Bohr";
        read_sync_double(input.relax_bfgs_init);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_axes");
        item.annotation = "which axes are fixed";
        read_sync_string(input.fixed_axes);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if ((para.input.fixed_axes == "shape" || para.input.fixed_axes == "volume") && !para.input.relax_new)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed shape and fixed volume only supported for relax_new = 1");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("fixed_ibrav");
        item.annotation = "whether to preseve lattice type during relaxation";
        read_sync_bool(input.fixed_ibrav);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.fixed_ibrav && !para.input.relax_new)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed_ibrav only available for relax_new = 1");
            }
            if (para.input.latname == "none" && para.input.fixed_ibrav)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "to use fixed_ibrav, latname must be provided");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("fixed_atoms");
        item.annotation = "whether to preseve direct coordinates of atoms "
                          "during relaxation";
        read_sync_bool(input.fixed_atoms);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.fixed_atoms && para.input.calculation == "relax")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "fixed_atoms is not meant to be used for calculation = relax");
            }
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO