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
                           "must still be commensurate with the ground-state k mesh. When "
                           "dfpt_compute_q0 is true, Gamma must be the FIRST q point and must "
                           "appear exactly once (the run rejects [X, Gamma] and [Gamma, X, "
                           "Gamma] orderings at input time, regardless of dfpt_out_dyn); "
                           "arbitrary orderings are allowed when dfpt_compute_q0 is false.";
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
                           "within the same DFPT run. Requires Gamma as the FIRST "
                           "q point, exactly once (the default dfpt_qmesh 1 1 1 "
                           "satisfies this); orderings such as [X, Gamma] or "
                           "[Gamma, X, Gamma] are rejected at input time, "
                           "independently of dfpt_out_dyn.";
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
        Input_Item item("dfpt_out_dyn");
        item.annotation = "write the raw dynamical matrix in QE dyn format "
                          "(pre-ASR, pre-LO-TO, unreduced force constants)";
        item.category = category;
        item.type = "Boolean";
        item.description = "Write one QE-compatible dyn file per q point into the output "
                           "directory: OUT.<suffix>/dyn.q<i>. Each file contains the "
                           "unreduced (pre-ASR, pre-LO-TO) force-constant matrix Phi in "
                           "Ry/bohr^2, the same content and layout as Quantum ESPRESSO's "
                           "ph.x output (fildyn), so the files can be diffed element-by-element "
                           "against QE references or post-processed by external tools. "
                           "At Gamma with dfpt_compute_q0 the Born effective charges and "
                           "dielectric tensor are appended in the QE convention.";
        item.default_value = "false";
        read_sync_bool(input.dfpt_out_dyn);
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_conv_thr");
        item.annotation = "convergence threshold of the DFPT response SCF";
        item.category = category;
        item.type = "Real";
        item.description = "Set the convergence threshold of the self-consistent DFPT cycle: "
                           "the iteration stops when the relative residual of the mixed "
                           "variable drops below this value for every displacement. With "
                           "dfpt_mix_type = v_sc_plain or v_sc_broyden (the default) the "
                           "residual is the screening-potential residual ||v_out - v_in|| / "
                           "||v_out||; with rho_plain, rho_kerker or rho_broyden it is the "
                           "first-order density residual ||drho_out - drho_in|| / ||drho_out||. "
                           "The electric-field (LO-TO) solve always mixes the density, using "
                           "the density-space counterpart of the selected scheme (e.g. "
                           "rho_broyden when dfpt_mix_type = v_sc_broyden), so its residual "
                           "is always the density residual.";
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
        item.domain.set_rule("value>=1");
        item.annotation = "max number of the DFPT first-order density mixing iterations";
        item.category = category;
        item.type = "Integer";
        item.description = "Set the maximum number of self-consistent DFPT iterations for "
                           "each atomic displacement.";
        item.default_value = "100";
        read_sync_int(input.dfpt_max_iter);
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_mix_beta");
        item.annotation = "mixing coefficient of the DFPT first-order density";
        item.category = category;
        item.type = "Real";
        item.description = "Set the plain-mixing coefficient of the DFPT response SCF. It "
                           "scales the preliminary step of the chosen dfpt_mix_type (the "
                           "plain mix itself, or the initial guess pushed into the Broyden "
                           "history). The response Jacobian has strongly negative eigenvalues "
                           "on the smallest-G shells (Coulomb stiffness), so beta must stay "
                           "below 2 / (1 + |lambda_min|); the default 0.4 keeps margin up to "
                           "|lambda_min| ~ 3. A larger value accelerates convergence for "
                           "weakly screened systems but may diverge.";
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
    {
        Input_Item item("dfpt_asr");
        item.annotation = "acoustic sum rule on the Gamma dynamical matrix: no, simple or crystal";
        item.category = category;
        item.type = "String";
        item.description = "Enforce the translational acoustic sum rule (sum over atoms of the "
                           "Gamma force-constant rows = 0) on the RAW Gamma dynamical matrix, "
                           "before the LO-TO nonanalytic term, together with the Z* sum rule "
                           "sum_a Z*_a = 0 on the Born effective charges. Unconstrained "
                           "linear-response self-consistency violates both at the percent to "
                           "tens-of-percent level (QE ph.x shows the same), and the reference "
                           "workflow post-processes with asr='crystal'. simple: subtract each "
                           "row sum from its diagonal element. crystal: the optimized "
                           "minimal-norm symmetric projection (keeps Hermiticity and "
                           "symmetry-equivalent atoms equivalent). no: leave the matrix as "
                           "computed.";
        item.default_value = "no";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            if (item.get_size() == 0)
            {
                para.input.dfpt_asr = "no";
            }
            else
            {
                para.input.dfpt_asr = strvalue;
            }
        };
        sync_string(input.dfpt_asr);
        item.check_value = [](const Input_Item& item, const Parameter& para) {
            if (para.input.dfpt_asr != "no" && para.input.dfpt_asr != "simple" && para.input.dfpt_asr != "crystal")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "dfpt_asr must be no, simple or crystal.");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("dfpt_mix_type");
        item.annotation = "mixing scheme of the DFPT response SCF: rho_plain, rho_kerker, rho_broyden, v_sc_plain or v_sc_broyden";
        item.category = category;
        item.type = "String";
        item.description = "Set the mixing scheme of the DFPT response SCF. rho_plain: the "
                           "first-order density is mixed with the same coefficient "
                           "dfpt_mix_beta on every shell. rho_kerker: the Kerker screen "
                           "f_g = |G+q|^2 / (|G+q|^2 + a^2) damps the Coulomb-stiff "
                           "smallest shells, stabilizing polar / small-gap systems at "
                           "beta up to 1; a^2 is set by dfpt_kerker_a2. rho_broyden: "
                           "Modified Broyden mixing (Johnson 1988) on the density, history "
                           "depth dfpt_broyden_ndim. v_sc_plain / v_sc_broyden: mix the "
                           "screening potential v_sc = v_H[drho] + v_xc[drho] instead of "
                           "the density, with plain or Modified Broyden coefficients; "
                           "v_sc_broyden is the default (the metric the Broyden least "
                           "squares works in is K-weighted, which stabilizes the "
                           "Coulomb-stiff shells that destabilize plain density mixing).";
         item.default_value = "v_sc_broyden";
         item.domain.set_choices({"rho_plain", "rho_kerker", "rho_broyden", "v_sc_plain", "v_sc_broyden"});
         // Contract: the line may be omitted (default_value v_sc_broyden
         // applies) or supply exactly one of the five choices. An
         // explicitly empty or unknown value is rejected by the scalar
         // choices domain before read_value runs; the previous
         // get_size()==0 fallback silently mapped an empty value onto
         // v_sc_broyden and is intentionally dropped.
         read_sync_string(input.dfpt_mix_type);
         this->add_item(item);
    }
    {
        Input_Item item("dfpt_kerker_a2");
        item.annotation = "Kerker screening parameter a^2 in 1/lat0^2 (dfpt_mix_type rho_kerker only)";
        item.category = category;
        item.type = "Real";
        item.description = "Set the Kerker screening parameter a^2 (in 1/lat0^2 units) used "
                           "by dfpt_mix_type = rho_kerker. The screen f_g = |G+q|^2 / "
                           "(|G+q|^2 + a^2) freezes the stiff long-wavelength shells; "
                           "larger a^2 damps more shells. Ignored by the other schemes.";
         item.default_value = "1.0";
         item.domain.set_rule("value>0");
         read_sync_double(input.dfpt_kerker_a2);
         this->add_item(item);
    }
    {
        Input_Item item("dfpt_broyden_ndim");
        item.annotation = "Broyden history depth for dfpt_mix_type = rho_broyden or v_sc_broyden";
        item.category = category;
        item.type = "Integer";
        item.description = "Set the number of history steps used by the "
                           "Modified Broyden mixing (Johnson 1988) when "
                           "dfpt_mix_type = rho_broyden or v_sc_broyden. Under "
                           "v_sc_broyden it sizes both Broyden instances: the "
                           "screening-potential mixer of the displacement loop "
                           "and the density mixer of the electric-field "
                           "(LO-TO) loop. More history generally accelerates "
                           "convergence at the cost of memory.";
         item.default_value = "8";
         item.domain.set_rule("value>=1");
         read_sync_int(input.dfpt_broyden_ndim);
         this->add_item(item);
    }
}
} // namespace ModuleIO
