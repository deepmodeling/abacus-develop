#include "read_input.h"
#include "read_input_tool.h"


namespace ModuleIO
{
void ReadInput::item_relax()
{
    {
        Input_Item item("ks_solver");
        item.annotation = "cg; dav; lapack; genelpa; scalapack_gvx; cusolver";
        read_sync_string(ks_solver);
        this->add_item(item);
    }
    {
        Input_Item item("scf_nmax");
        item.annotation = "number of electron iterations";
        read_sync_int(scf_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_nmax");
        item.annotation = "number of ion iteration steps";
        read_sync_int(relax_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("out_stru");
        item.annotation = "output the structure files after each ion step";
        read_sync_bool(out_stru);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr");
        item.annotation = "force threshold, unit: Ry/Bohr";
        read_sync_double(force_thr);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev");
        item.annotation = "force threshold, unit: eV/Angstrom";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.force_thr_ev = doublevalue ;
        };
        item.resetvalue = [](const Input_Item& item, Parameter& para) {
            para.force_thr = para.force_thr_ev / 13.6058 * 0.529177;
        };
        sync_double(force_thr_ev);
        this->add_item(item);
    }
    {
        Input_Item item("force_thr_ev2");
        item.annotation = "force invalid threshold, unit: eV/Angstrom";
        read_sync_double(force_thr_ev2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_cg_thr");
        item.annotation = "threshold for switching from cg to bfgs, unit: eV/Angstrom";
        read_sync_double(relax_cg_thr);
        this->add_item(item);
    }
    {
        Input_Item item("stress_thr");
        item.annotation = "stress threshold";
        read_sync_double(stress_thr);
        this->add_item(item);
    }
    {
        Input_Item item("press1");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(press1);
        this->add_item(item);
    }
    {
        Input_Item item("press2");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(press2);
        this->add_item(item);
    }
    {
        Input_Item item("press3");
        item.annotation = "target pressure, unit: KBar";
        read_sync_double(press3);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w1");
        item.annotation = "wolfe condition 1 for bfgs";
        read_sync_double(relax_bfgs_w1);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_w2");
        item.annotation = "wolfe condition 2 for bfgs";
        read_sync_double(relax_bfgs_w2);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmax");
        item.annotation = "maximal trust radius, unit: Bohr";
        read_sync_double(relax_bfgs_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_rmin");
        item.annotation = "minimal trust radius, unit: Bohr";
        read_sync_double(relax_bfgs_rmin);
        this->add_item(item);
    }
    {
        Input_Item item("relax_bfgs_init");
        item.annotation = "initial trust radius, unit: Bohr";
        read_sync_double(relax_bfgs_init);
        this->add_item(item);
    }
    {
        Input_Item item("cal_stress");
        item.annotation = "calculate the stress or not";
        read_sync_bool(cal_stress);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_axes");
        item.annotation = "which axes are fixed";
        read_sync_string(fixed_axes);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_ibrav");
        item.annotation = "whether to preseve lattice type during relaxation";
        read_sync_bool(fixed_ibrav);
        this->add_item(item);
    }
    {
        Input_Item item("fixed_atoms");
        item.annotation = "whether to preseve direct coordinates of atoms during relaxation";
        read_sync_bool(fixed_atoms);
        this->add_item(item);
    }
    {
        Input_Item item("relax_method");
        item.annotation = "cg; bfgs; sd; cg; cg_bfgs;";
        read_sync_string(relax_method);
        this->add_item(item);
    }
    {
        Input_Item item("relax_new");
        item.annotation = "whether to use the new relaxation method";
        read_sync_bool(relax_new);
        this->add_item(item);
    }
    {
        Input_Item item("relax_scale_force");
        item.annotation = "controls the size of the first CG step if relax_new is true";
        read_sync_double(relax_scale_force);
        this->add_item(item);
    }
    {
        Input_Item item("out_level");
        item.annotation = "ie(for electrons); i(for ions);";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.out_level = strvalue;
            para.out_md_control = true;
        };
        sync_string(out_level);
        add_bool_bcast(out_md_control); //Since "out_md_control" has been assigned a value, it needs to be broadcasted
        this->add_item(item);
    }
    {
        Input_Item item("out_dm");
        item.annotation = ">0 output density matrix";
        read_sync_bool(out_dm);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm1");
        item.annotation = ">0 output density matrix (multi-k points)";
        read_sync_bool(out_dm1);
        this->add_item(item);
    }
    {
        Input_Item item("out_bandgap");
        item.annotation = "if true, print out bandgap";
        read_sync_bool(out_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("use_paw");
        item.annotation = "whether to use PAW in pw calculation";
        read_sync_bool(use_paw);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_labels");
        item.annotation = ">0 compute descriptor for deepks";
        read_sync_bool(deepks_out_labels);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_scf");
        item.annotation = ">0 add V_delta to Hamiltonian";
        read_sync_bool(deepks_scf);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_equiv");
        item.annotation = "whether to use equivariant version of DeePKS";
        read_sync_bool(deepks_equiv);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_bandgap");
        item.annotation = ">0 for bandgap label";
        read_sync_bool(deepks_bandgap);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_out_unittest");
        item.annotation = "if set 1, prints intermediate quantities that shall be used for making unit test";
        read_sync_bool(deepks_out_unittest);
        this->add_item(item);
    }
    {
        Input_Item item("deepks_model");
        item.annotation = "file dir of traced pytorch model: 'model.ptg";
        read_sync_string(deepks_model);
        this->add_item(item);
    }
}
} // namespace ModuleIO