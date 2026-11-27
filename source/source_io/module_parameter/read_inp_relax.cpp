#include "source_base/global_function.h"
#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

namespace ModuleIO
{

void ReadInput::item_relax()
{
    // NOTE: The order of add_item() calls below determines the parameter order
    // in the generated documentation (docs/advanced/input_files/input-main.md).
    // Please preserve this ordering when adding new parameters.
    {
        Input_Item item("relax_method");
        item.domain.set_choices({"cg", "cg_trad", "bfgs", "bfgs_trad", "sd", "cg_bfgs", "lbfgs"});
        item.annotation = "cg; cg_trad; bfgs; bfgs_trad; sd; cg_bfgs; lbfgs";
        item.category = "Geometry relaxation";
        item.type = "String";
        item.description = R"(The method used for geometry optimization. Select exactly one option:

* cg: Simultaneously optimize ionic positions and cell parameters with conjugate gradient and line search (recommended default).
* cg_trad: First optimize ionic positions at fixed cell, then update the cell, and repeat.
* bfgs: Default BFGS that directly updates the inverse Hessian (recommended).
* bfgs_trad: Traditional BFGS that updates the Hessian matrix B and then inverts it.
* lbfgs: Limited-memory BFGS algorithm, suitable for large systems.
* cg_bfgs: Mixed method starting with CG and switching to BFGS when force convergence reaches relax_cg_thr.
* sd: Steepest descent algorithm. Not recommended for production use.)";
        item.default_value = "cg";
        item.unit = "";
        read_sync_string(input.relax_method);
        this->add_item(item);
    }
    {
        Input_Item item("relax_scale_force");
        item.annotation = "controls the size of the first simultaneous CG step";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The paramether controls the size of the first conjugate gradient step. A smaller value means the first step along a new CG direction is smaller. This might be helpful for large systems, where it is safer to take a smaller initial step to prevent the collapse of the whole configuration.";
        item.default_value = "0.5";
        item.unit = "";
        item.set_availability("relax_method==cg");
        read_sync_double(input.relax_scale_force);
        this->add_item(item);
    }
    {
        Input_Item item("relax_nmax");
        item.annotation = "number of ion iteration steps";
        item.category = "Geometry relaxation";
        item.type = "Integer";
        item.description = "The maximal number of ionic iteration steps. If set to 0, the code performs a quick \"dry run\", stopping just after initialization. This is useful to check for input correctness and to have the summary printed.";
        item.default_value = "1 for SCF, 50 for relax and cell-relax calcualtions";
        item.unit = "";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            const std::string& calculation = para.input.calculation;
            const std::vector<std::string> singlelist
                = {"scf", "nscf", "get_s", "get_pchg", "get_wf", "test_memory", "test_neighbour", "gen_bessel"};
            if (std::find(singlelist.begin(), singlelist.end(), calculation) != singlelist.end())
            {
                if (para.input.relax_nmax != 0)
                {
                    para.input.relax_nmax = 1;
                }
            }
            else if (calculation == "relax" || calculation == "cell-relax")
            {
                if (para.input.relax_nmax < 0)
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
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "When relax_method is set to cg_bfgs, a mixed algorithm of conjugate gradient (CG) and Broyden–Fletcher–Goldfarb–Shanno (BFGS) is used. The ions first move according to the CG method, then switch to the BFGS method when the maximum force on atoms is reduced below this threshold.";
        item.default_value = "0.5";
        item.unit = "eV/Angstrom";
        item.set_availability("relax_method==cg_bfgs");
        read_sync_double(input.relax_cg_thr);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr");
        item.domain.set_rule("value>=0");
        item.annotation = "force threshold, unit: Ry/Bohr";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Threshold of the force convergence. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to force_thr_ev except for the unit, you can choose either you like.";
        item.default_value = "0.001";
        item.unit = "Ry/Bohr (25.7112 eV/Angstrom)";
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
        item.domain.set_rule("value>=0");
        item.annotation = "force threshold, unit: eV/Angstrom";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Threshold of the force convergence. The threshold is compared with the largest force among all of the atoms. The recommended value for using atomic orbitals is 0.04 eV/Angstrom (0.0016 Ry/Bohr). The parameter is equivalent to force_thr except for the unit. You may choose either you like.";
        item.default_value = "0.0257112";
        item.unit = "eV/Angstrom (0.03889 Ry/Bohr)";
        item.read_value = [](const Input_Item& item, Parameter& para) { para.input.force_thr_ev = doublevalue; };
        sync_double(input.force_thr_ev);
        this->add_item(item);
    }
    {
        Input_Item item("force_zero_out");
        item.annotation = "force invalid threshold, unit: eV/Angstrom";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The atomic forces that are smaller than force_zero_out will be treated as zero.";
        item.default_value = "0.0";
        item.unit = "eV/Angstrom";
        read_sync_double(input.force_zero_out);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w1");
        item.annotation = "wolfe condition 1 for bfgs";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Controls the Wolfe condition for the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. This parameter sets the sufficient decrease condition (c1 in Wolfe conditions). For more information, see Phys. Chem. Chem. Phys., 2000, 2, 2177.";
        item.default_value = "0.01";
        item.unit = "";
        item.set_availability("relax_method in [bfgs, bfgs_trad, cg_bfgs]");
        read_sync_double(input.relax_bfgs_w1);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w2");
        item.annotation = "wolfe condition 2 for bfgs";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Controls the Wolfe condition for the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm used in geometry relaxation. This parameter sets the curvature condition (c2 in Wolfe conditions). For more information, see Phys. Chem. Chem. Phys., 2000, 2, 2177.";
        item.default_value = "0.5";
        item.unit = "";
        item.set_availability("relax_method in [bfgs, bfgs_trad, cg_bfgs]");
        read_sync_double(input.relax_bfgs_w2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmax");
        item.annotation = "maximal trust radius, unit: Bohr";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Maximum allowed total displacement of all atoms during geometry optimization. The sum of atomic displacements can increase during optimization steps but cannot exceed this value.";
        item.default_value = "0.8";
        item.unit = "Bohr";
        item.set_availability("relax_method in [bfgs, bfgs_trad, cg_bfgs]");
        read_sync_double(input.relax_bfgs_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmin");
        item.annotation = "minimal trust radius, unit: Bohr";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Minimum allowed total displacement of all atoms. When the total atomic displacement falls below this value and force convergence is not achieved, the calculation will terminate. Note: This parameter is not used in the default BFGS algorithm (relax_method = bfgs).";
        item.default_value = "1e-5";
        item.unit = "Bohr";
        item.set_availability("relax_method==bfgs_trad");
        read_sync_double(input.relax_bfgs_rmin);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_init");
        item.annotation = "initial trust radius, unit: Bohr";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "Initial total displacement of all atoms in the first BFGS step. This sets the scale for the initial movement.";
        item.default_value = "0.5";
        item.unit = "Bohr";
        item.set_availability("relax_method in [bfgs, bfgs_trad, cg_bfgs]");
        read_sync_double(input.relax_bfgs_init);
        this->add_item(item);
    }
    {
        Input_Item item("stress_thr");
        item.annotation = "stress threshold";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The threshold of the stress convergence. The threshold is compared with the largest component of the stress tensor.";
        item.default_value = "0.5";
        item.unit = "kbar";
        read_sync_double(input.stress_thr);
        this->add_item(item);
    }
    {
        Input_Item item("press1");
        item.annotation = "target pressure, unit: KBar";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The external pressures along three axes. Positive input value is taken as compressive stress.";
        item.default_value = "0";
        item.unit = "kbar";
        read_sync_double(input.press1);
        this->add_item(item);
    }
    {
        Input_Item item("press2");
        item.annotation = "target pressure, unit: KBar";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The external pressures along three axes. Positive input value is taken as compressive stress.";
        item.default_value = "0";
        item.unit = "kbar";
        read_sync_double(input.press2);
        this->add_item(item);
    }
    {
        Input_Item item("press3");
        item.annotation = "target pressure, unit: KBar";
        item.category = "Geometry relaxation";
        item.type = "Real";
        item.description = "The external pressures along three axes. Positive input value is taken as compressive stress.";
        item.default_value = "0";
        item.unit = "kbar";
        read_sync_double(input.press3);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_axes");
        item.domain.add_requirement("fixed_axes in [shape, volume]", "relax_method==cg");
        item.annotation = "which axes are fixed";
        item.category = "Geometry relaxation";
        item.type = "String";
        item.description = R"(Specifies which cell degrees of freedom are fixed during variable-cell relaxation. The available options depend on relax_method:

With relax_method = cg (default), all options are available:
* None: Default; all cell parameters can relax freely
* volume: Relaxation with fixed volume (allows shape changes)
* shape: Fix shape but allow volume changes (hydrostatic pressure only)
* a: Fix the a-axis lattice vector during relaxation
* b: Fix the b-axis lattice vector during relaxation
* c: Fix the c-axis lattice vector during relaxation
* ab: Fix both a and b axes during relaxation
* ac: Fix both a and c axes during relaxation
* bc: Fix both b and c axes during relaxation
* abc: Fix all three lattice vectors during relaxation

With relax_method set to cg_trad, bfgs, bfgs_trad, lbfgs, sd, or cg_bfgs, None and a, b, c, ab, ac, bc, abc are available. The shape and volume options require cg.

[NOTE] For VASP users, see the ISIF correspondence table in the geometry optimization documentation.)";
        item.default_value = "None";
        item.unit = "";
        item.set_availability("calculation==cell-relax");
        read_sync_string(input.fixed_axes);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_ibrav");
        item.domain.add_requirement("fixed_ibrav==true", "relax_method==cg and latname!=none");
        item.annotation = "whether to preseve lattice type during relaxation";
        item.category = "Geometry relaxation";
        item.type = "Boolean";
        item.description = R"(* True: the lattice type will be preserved during relaxation. The lattice vectors are reconstructed to match the specified Bravais lattice type after each update.
* False: No restrictions are exerted during relaxation in terms of lattice type

[NOTE] Note: it is possible to use fixed_ibrav with fixed_axes, but please make sure you know what you are doing. For example, if we are doing relaxation of a simple cubic lattice (latname = "sc"), and we use fixed_ibrav along with fixed_axes = "volume", then the cell is never allowed to move and as a result, the relaxation never converges. When both are used, fixed_ibrav is applied first, then fixed_axes = "volume" rescaling is applied.)";
        item.default_value = "False";
        item.unit = "";
        item.set_availability("relax_method==cg and latname!=none");
        read_sync_bool(input.fixed_ibrav);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_atoms");
        item.domain.add_requirement("fixed_atoms==true", "calculation!=relax");
        item.annotation = "whether to preseve direct coordinates of atoms "
                          "during relaxation";
        item.category = "Geometry relaxation";
        item.type = "Boolean";
        item.description = R"(* True: The direct coordinates of atoms will be preserved during variable-cell relaxation.
* False: No restrictions are exerted on positions of all atoms. However, users can still fix certain components of certain atoms by using the m keyword in STRU file. For the latter option, check the end of this instruction.)";
        item.default_value = "False";
        item.unit = "";
        read_sync_bool(input.fixed_atoms);
        this->add_item(item);
    }
    {
        Input_Item item("relax");
        item.annotation = "allow relaxation along the specific direction";
        read_sync_bool(input.relax);
        this->add_item(item);
    }
}
} // namespace ModuleIO
