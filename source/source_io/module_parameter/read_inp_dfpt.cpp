#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>

namespace ModuleIO
{
void ReadInput::item_dfpt()
{
    const std::string category = "Density functional perturbation theory";
    {
        Input_Item item("dfpt_qmesh");
        item.annotation = "Monkhorst-Pack q mesh for DFPT, should be >= 1, "
                          "and commensurate with the KPT mesh";
        item.category = category;
        item.type = "Vector of Int (1 or 3 values)";
        item.description = "Set the Monkhorst-Pack q mesh (gamma-centered) for DFPT phonon "
                           "calculations. The q mesh must be commensurate with the ground-state "
                           "k mesh: k + q must be a point of the k list (modulo a reciprocal "
                           "lattice vector). For example, a 4x4x4 KPT mesh is commensurate with "
                           "dfpt_qmesh values of 1, 2, or 4 along each direction. This parameter "
                           "is ignored when dfpt_qfile is set.";
        item.default_value = "1 1 1";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            if (count == 1)
            {
                para.input.dfpt_qmesh[0] = para.input.dfpt_qmesh[1] = para.input.dfpt_qmesh[2] = intvalue;
            }
            else if (count == 3)
            {
                para.input.dfpt_qmesh[0] = std::stoi(item.str_values[0]);
                para.input.dfpt_qmesh[1] = std::stoi(item.str_values[1]);
                para.input.dfpt_qmesh[2] = std::stoi(item.str_values[2]);
            }
            else
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_qmesh can only accept one or three values.");
            }
        };
        sync_intvec(input.dfpt_qmesh, 3, 1);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            for (int i = 0; i < 3; i++)
            {
                if (para.input.dfpt_qmesh[i] < 1)
                {
                    ModuleBase::WARNING_QUIT("ReadInput", "dfpt_qmesh must be >= 1.");
                }
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_qfile");
        item.annotation = "file containing the DFPT q points; empty means the "
                          "dfpt_qmesh Monkhorst-Pack mesh";
        item.category = category;
        item.type = "String";
        item.description = "Set the file containing the q points for DFPT, in the same format "
                           "as the KPT file (Q_POINTS card: Gamma/Monkhorst-Pack mesh, or an "
                           "explicit Direct/Cartesian list; symmetry reduction is not applied "
                           "to file q lists). When set, it overrides dfpt_qmesh. Each q point "
                           "must still be commensurate with the ground-state k mesh.";
        item.default_value = "\"\"";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            if (item.get_size() == 0)
            {
                para.input.dfpt_qfile = "";
            }
            else
            {
                para.input.dfpt_qfile = strvalue;
            }
        };
        sync_string(input.dfpt_qfile);
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_compute_q0");
        item.annotation = "whether to compute epsilon_inf and Born effective "
                          "charges at q = 0";
        item.category = category;
        item.type = "Boolean";
        item.description = "Whether to compute the macroscopic dielectric tensor "
                           "(epsilon_inf) and the Born effective charges at q = 0 "
                           "within the same DFPT run. Requires a q point at Gamma "
                           "(the default dfpt_qmesh 1 1 1).";
        item.default_value = "false";
        read_sync_bool(input.dfpt_compute_q0);
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_loto");
        item.annotation = "whether to apply the LO-TO non-analytic correction "
                          "at q = 0";
        item.category = category;
        item.type = "Boolean";
        item.description = "Whether to apply the Lyddane-Sachs-Teller non-analytic correction "
                           "to the Gamma-point dynamical matrix, which splits the longitudinal "
                           "and transverse optical modes. Requires dfpt_compute_q0 to be true, "
                           "since the correction is built from epsilon_inf and the Born "
                           "effective charges.";
        item.default_value = "false";
        read_sync_bool(input.dfpt_loto);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dfpt_loto && !para.input.dfpt_compute_q0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_loto requires dfpt_compute_q0 to be true.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_conv_thr");
        item.annotation = "convergence threshold of the DFPT first-order density";
        item.category = category;
        item.type = "Real";
        item.description = "Set the convergence threshold of the self-consistent DFPT cycle: "
                           "the iteration stops when the relative residual of the first-order "
                           "density ||drho_out - drho_in|| / ||drho_out|| drops below this "
                           "value for every displacement.";
        item.default_value = "1.0e-8";
        read_sync_double(input.dfpt_conv_thr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dfpt_conv_thr <= 0.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_conv_thr must be > 0.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_max_iter");
        item.annotation = "max number of the DFPT first-order density mixing iterations";
        item.category = category;
        item.type = "Integer";
        item.description = "Set the maximum number of self-consistent DFPT iterations for "
                           "each atomic displacement.";
        item.default_value = "100";
        read_sync_int(input.dfpt_max_iter);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dfpt_max_iter < 1)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_max_iter must be >= 1.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_mix_beta");
        item.annotation = "mixing coefficient of the DFPT first-order density";
        item.category = category;
        item.type = "Real";
        item.description = "Set the plain-mixing coefficient of the first-order density in the "
                           "self-consistent DFPT cycle. The response Jacobian has strongly "
                           "negative eigenvalues on the smallest-G shells (Coulomb stiffness), "
                           "so beta must stay below 2 / (1 + |lambda_min|); the default 0.4 "
                           "keeps margin up to |lambda_min| ~ 3. A larger value accelerates "
                           "convergence for weakly screened systems but may diverge.";
        item.default_value = "0.4";
        read_sync_double(input.dfpt_mix_beta);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dfpt_mix_beta <= 0.0 || para.input.dfpt_mix_beta > 1.0)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_mix_beta must be in (0, 1].");
            }
        };
        this->add_item(item);
    }
}
} // namespace ModuleIO
