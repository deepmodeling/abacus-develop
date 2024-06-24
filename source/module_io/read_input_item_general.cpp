#include "read_input.h"
#include "read_input_tool.h"

#include <fstream>

namespace ModuleIO
{
// There are some examples:
// Generallly:
// {
//      Input_Item item("suffix");
//      item.annotation = "the name of main output directory";
//      read_sync_string(suffix);
//      this->add_item(item);
// }
//
// Specially:
// {
//      Input_Item item("kspacing");
//      item.annotation = "unit in 1/bohr, should be > 0, default is 0 which means read KPT file";
//
//      item.readvalue = [](const Input_Item& item, Parameter& para) {
//          para.input.kspacing[0] = convertstr<double>(item.str_values[0]);
//          para.input.kspacing[1] = convertstr<double>(item.str_values[1]);
//          para.input.kspacing[2] = convertstr<double>(item.str_values[2]);
//      };
//
//      item.resetvalue = [](const Input_Item& item, Parameter& para) {
//          if(para.input.kspacing[0] <= 0) para.input.kspacing[0] = 1;
//      };
//
//      item.checkvalue = [](const Input_Item& item, const Parameter& para) {assert(para.input.kspacing[0]>0);};
//
//      item.getfinalvalue = [](Input_Item& item, const Parameter& para) {
//          item.final_value << para.input.kspacing[0] << " " << para.input.kspacing[1] << " " <<
//          para.input.kspacing[2];
//      };
//
//      add_doublevec_bcast(&Parameter::PARAMETER, N);
//      this->add_item(item);
//  }
void ReadInput::item_general()
{
    {
        Input_Item item("suffix");
        item.annotation = "the name of main output directory";
        read_sync_string(suffix);
        this->add_item(item);
    }
    {
        Input_Item item("latname");
        item.annotation = "the name of lattice name";
        read_sync_string(latname);
        this->add_item(item);
    }
    {
        Input_Item item("stru_file");
        item.annotation = "the filename of file containing atom positions";
        read_sync_string(stru_file);
        this->add_item(item);
    }
    {
        Input_Item item("kpoint_file");
        item.annotation = "the name of file containing k points";
        read_sync_string(kpoint_file);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_dir");
        item.annotation = "the directory containing pseudo files";
        read_sync_string(pseudo_dir);
        this->add_item(item);
    }
    {
        Input_Item item("orbital_dir");
        item.annotation = "the directory containing orbital files";
        read_sync_string(orbital_dir);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_rcut");
        item.annotation = "default #exchange correlation functional";
        read_sync_double(pseudo_rcut);
        this->add_item(item);
    }
    {
        Input_Item item("pseudo_mesh");
        item.annotation = "0: use our own mesh to do radial renormalization; 1: use mesh as in QE";
        read_sync_bool(pseudo_mesh);
        this->add_item(item);
    }
    {
        Input_Item item("lmaxmax");
        item.annotation = "maximum of l channels used";
        read_sync_int(lmaxmax);
        this->add_item(item);
    }
    {
        Input_Item item("dft_functional");
        item.annotation = "exchange correlation functional";
        read_sync_string(dft_functional);
        this->add_item(item);
    }
    {
        Input_Item item("xc_temperature");
        item.annotation = "temperature for finite temperature functionals";
        read_sync_double(xc_temperature);
        this->add_item(item);
    }
    {
        Input_Item item("calculation");
        item.annotation = "test; scf; relax; nscf; get_wf; get_pchg";
        read_sync_string(calculation);
        this->add_item(item);
    }
    {
        Input_Item item("esolver_type");
        item.annotation = "the energy solver: ksdft, sdft, ofdft, tddft, lj, dp";
        read_sync_string(esolver_type);
        this->add_item(item);
    }
    {
        Input_Item item("ntype");
        item.annotation = "atom species number";
        read_sync_int(ntype);
        this->add_item(item);
    }
    {
        Input_Item item("nspin");
        item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
        read_sync_int(nspin);
        this->add_item(item);
    }
    {
        Input_Item item("kspacing");
        item.annotation = "unit in 1/bohr, should be > 0, default is 0 which means read KPT file";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            if (count == 1)
            {
                para.input.kspacing[0] = para.input.kspacing[1] = para.input.kspacing[2] = doublevalue;
            }
            else if (count == 3)
            {
                para.input.kspacing[0] = convertstr<double>(item.str_values[0]);
                para.input.kspacing[1] = convertstr<double>(item.str_values[1]);
                para.input.kspacing[2] = convertstr<double>(item.str_values[2]);
            }
            else
            {
                throw std::runtime_error("kspacing can only accept one or three double values.");
            }
        };
        sync_doublevec(kspacing, 3);
        this->add_item(item);
    }
    {
        Input_Item item("min_dist_coef");
        item.annotation = "factor related to the allowed minimum distance between two atoms";
        read_sync_double(min_dist_coef);
        this->add_item(item);
    }
    {
        Input_Item item("nbands");
        item.annotation = "number of bands";
        read_sync_int(nbands);
        this->add_item(item);
    }
    {
        Input_Item item("nbands_istate");
        item.annotation = "number of bands around Fermi level for get_pchg calulation";
        read_sync_int(nbands_istate);
        this->add_item(item);
    }
    {
        Input_Item item("bands_to_print");
        item.annotation = "specify the bands to be calculated in the get_pchg calculation";
        read_sync_string(bands_to_print);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry");
        item.annotation = "the control of symmetry";
        read_sync_string(symmetry);
        this->add_item(item);
    }
    {
        Input_Item item("init_vel");
        item.annotation = "read velocity from STRU or not";
        read_sync_bool(init_vel);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_prec");
        item.annotation = "accuracy for symmetry";
        read_sync_double(symmetry_prec);
        this->add_item(item);
    }
    {
        Input_Item item("symmetry_autoclose");
        item.annotation = "whether to close symmetry automatically when error occurs in symmetry analysis";
        read_sync_bool(symmetry_autoclose);
        this->add_item(item);
    }
    {
        Input_Item item("nelec");
        item.annotation = "input number of electrons";
        read_sync_double(nelec);
        this->add_item(item);
    }
    {
        Input_Item item("nelec_delta");
        item.annotation = "change in the number of total electrons";
        read_sync_double(nelec_delta);
        this->add_item(item);
    }
    {
        Input_Item item("nupdown");
        item.annotation = "the difference number of electrons between spin-up and spin-down";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.input.nupdown = doublevalue;
            para.input.two_fermi = true;
        };

        sync_double(nupdown);
        add_bool_bcast(two_fermi); // two_fermi is also changed, so need to bcast
        this->add_item(item);
    }
    {
        Input_Item item("out_mul");
        item.annotation = "mulliken charge or not";
        read_sync_bool(out_mul);
        this->add_item(item);
    }
    {
        Input_Item item("noncolin");
        item.annotation = "using non-collinear-spin";
        read_sync_bool(noncolin);
        this->add_item(item);
    }
    {
        Input_Item item("lspinorb");
        item.annotation = "consider the spin-orbit interaction";
        read_sync_bool(lspinorb);
        this->add_item(item);
    }
    {
        Input_Item item("kpar");
        item.annotation = "devide all processors into kpar groups and k points will be distributed among";
        read_sync_int(kpar);
        this->add_item(item);
    }
    {
        Input_Item item("bndpar");
        item.annotation = "devide all processors into bndpar groups and bands will be distributed among each group";
        read_sync_int(bndpar);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_elec");
        item.annotation = "the frequency ( >= 0) of electronic iter to output charge density and wavefunction. 0: "
                          "output only when converged";
        read_sync_int(out_freq_elec);
        this->add_item(item);
    }
    {
        Input_Item item("dft_plus_dmft");
        item.annotation = "true:DFT+DMFT; false: standard DFT calcullation(default)";
        read_sync_bool(dft_plus_dmft);
        this->add_item(item);
    }
    {
        Input_Item item("rpa");
        item.annotation = "true:generate output files used in rpa calculation; false:(default)";
        read_sync_bool(rpa);
        this->add_item(item);
    }
    {
        Input_Item item("printe");
        item.annotation = "Print out energy for each band for every printe steps";
        read_sync_int(printe);
        this->add_item(item);
    }
    {
        Input_Item item("mem_saver");
        item.annotation = "Only for nscf calculations. if set to 1, then a memory saving technique will be used for "
                          "many k point calculations.";
        read_sync_int(mem_saver);
        this->add_item(item);
    }
    {
        Input_Item item("diago_proc");
        item.annotation = "the number of procs used to do diagonalization";
        read_sync_int(diago_proc);
        this->add_item(item);
    }
    {
        Input_Item item("nbspline");
        item.annotation = "the order of B-spline basis";
        read_sync_int(nbspline);
        this->add_item(item);
    }
    {
        Input_Item item("wannier_card");
        item.annotation = "input card for wannier functions";
        read_sync_string(wannier_card);
        this->add_item(item);
    }
    {
        Input_Item item("soc_lambda");
        item.annotation = "The fraction of averaged SOC pseudopotential is given by (1-soc_lambda)";
        read_sync_double(soc_lambda);
        this->add_item(item);
    }
    {
        Input_Item item("cal_force");
        item.annotation = "if calculate the force at the end of the electronic iteration";
        read_sync_bool(cal_force);
        this->add_item(item);
    }
    {
        Input_Item item("out_freq_ion");
        item.annotation = "the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output "
                          "only when ion steps are finished";
        read_sync_int(out_freq_ion);
        this->add_item(item);
    }
    {
        Input_Item item("elpa_num_thread");
        item.annotation = "Number of threads need to use in elpa";
        read_sync_int(elpa_num_thread);
        this->add_item(item);
    }
    {
        Input_Item item("device");
        item.annotation = "the computing device for ABACUS";
        read_sync_string(device);
        this->add_item(item);
    }
    {
        Input_Item item("precision");
        item.annotation = "the computing precision for ABACUS";
        read_sync_string(precision);
        this->add_item(item);
    }
}

} // namespace ModuleIO