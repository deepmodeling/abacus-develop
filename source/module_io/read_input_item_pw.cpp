#include "read_input.h"
#include "read_input_tool.h"


namespace ModuleIO
{
void ReadInput::item_pw()
{
    {
        Input_Item item("ecutwfc");
        item.annotation = "energy cutoff for wave functions";
        read_bcast_double(ecutwfc);
        this->add_item(item);
    }
    {
        Input_Item item("ecutrho");
        item.annotation = "energy cutoff for charge density and potential";
        read_bcast_double(ecutrho);
        this->add_item(item);
    }
    {
        Input_Item item("erf_ecut");
        item.annotation = "the value of the constant energy cutoff";
        read_bcast_double(erf_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("erf_height");
        item.annotation = "the height of the energy step for reciprocal vectors";
        read_bcast_double(erf_height);
        this->add_item(item);
    }
    {
        Input_Item item("erf_sigma");
        item.annotation = "the width of the energy step for reciprocal vectors";
        read_bcast_double(erf_sigma);
        this->add_item(item);
    }
    {
        Input_Item item("fft_mode");
        item.annotation = "mode of FFTW";
        read_bcast_int(fft_mode);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_nmax");
        item.annotation = "max iteration number for cg";
        read_bcast_int(pw_diag_nmax);
        this->add_item(item);
    }
    {
        Input_Item item("diago_cg_prec");
        item.annotation = "diago_cg_prec";
        read_bcast_int(diago_cg_prec);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_ndim");
        item.annotation = "dimension of workspace for Davidson diagonalization";
        read_bcast_int(pw_diag_ndim);
        this->add_item(item);
    }
    {
        Input_Item item("diago_full_acc");
        item.annotation = "all the empty states are diagonalized";
        read_bcast_bool(diago_full_acc);
        this->add_item(item);
    }
    {
        Input_Item item("pw_diag_thr");
        item.annotation = "threshold for eigenvalues is cg electron iterations";
        read_bcast_double(pw_diag_thr);
        this->add_item(item);
    }
    {
        Input_Item item("nb2d");
        item.annotation = "matrix 2d division";
        read_bcast_int(nb2d);
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr");
        item.annotation = "charge density error";
        read_bcast_double(scf_thr);
        this->add_item(item);
    }
    {
        Input_Item item("scf_thr_type");
        item.annotation = "type of the criterion of scf_thr, 1: reci drho for pw, 2: real drho for lcao";
        read_bcast_int(scf_thr_type);
        this->add_item(item);
    }
    {
        Input_Item item("init_wfc");
        item.annotation = "start wave functions are from 'atomic', 'atomic+random', 'random' or";
        read_bcast_string(init_wfc);
        this->add_item(item);
    }
    {
        Input_Item item("psi_initializer");
        item.annotation = "whether to use psi_initializer";
        read_bcast_bool(psi_initializer);
        this->add_item(item);
    }
    {
        Input_Item item("init_chg");
        item.annotation = "start charge is from 'atomic' or file";
        read_bcast_string(init_chg);
        this->add_item(item);
    }
    {
        Input_Item item("chg_extrap");
        item.annotation = "atomic; first-order; second-order; dm:coefficients of SIA";
        read_bcast_string(chg_extrap);
        this->add_item(item);
    }
    {
        Input_Item item("out_chg");
        item.annotation = ">0 output charge density for selected electron steps";
        read_bcast_int(out_chg);
        this->add_item(item);
    }
    {
        Input_Item item("out_pot");
        item.annotation = "output realspace potential";
        read_bcast_int(out_pot);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_pw");
        item.annotation = "output wave functions";
        read_bcast_int(out_wfc_pw);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_r");
        item.annotation = "output wave functions in realspace";
        read_bcast_bool(out_wfc_r);
        this->add_item(item);
    }
    {
        Input_Item item("out_dos");
        item.annotation = "output energy and dos";
        read_bcast_int(out_dos);
        this->add_item(item);
    }
    {
        Input_Item item("out_band");
        item.annotation = "output energy and band structure (with precision 8)";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            if (count == 1)
            {
                para.out_band[0] = convertstr<int>(item.str_values[0]);
                para.out_band[1] = 8;
            }
            else if (count == 2)
            {
                para.out_band[0] = convertstr<int>(item.str_values[0]);
                para.out_band[1] = convertstr<int>(item.str_values[1]);
            }
            else
            {
                throw std::runtime_error("out_band should have 1 or 2 values");
            }
        };
        bcast_intvec_param(out_band.data(), 2);
        this->add_item(item);
    }
    {
        Input_Item item("out_proj_band");
        item.annotation = "output projected band structure";
        read_bcast_bool(out_proj_band);
        this->add_item(item);
    }
    {
        Input_Item item("restart_save");
        item.annotation = "print to disk every step for restart";
        read_bcast_bool(restart_save);
        this->add_item(item);
    }
    {
        Input_Item item("restart_load");
        item.annotation = "restart from disk";
        read_bcast_bool(restart_load);
        this->add_item(item);
    }
    {
        Input_Item item("read_file_dir");
        item.annotation = "directory of files for reading";
        read_bcast_string(read_file_dir);
        this->add_item(item);
    }
    {
        Input_Item item("nx");
        item.annotation = "number of points along x axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.nx = intvalue;
            para.ncx = intvalue;
        };
        bcast_int_param(nx);
        bcast_int_param(ncx);
        this->add_item(item);
    }
    {
        Input_Item item("ny");
        item.annotation = "number of points along y axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.ny = intvalue;
            para.ncy = intvalue;
        };
        bcast_int_param(ny);
        bcast_int_param(ncy);
        this->add_item(item);
    }
    {
        Input_Item item("nz");
        item.annotation = "number of points along z axis for FFT grid";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            para.nz = intvalue;
            para.ncz = intvalue;
        };
        bcast_int_param(nz);
        bcast_int_param(ncz);
        this->add_item(item);
    }
    {
        Input_Item item("ndx");
        item.annotation = "number of points along x axis for FFT smooth grid";
        read_bcast_int(ndx);
        this->add_item(item);
    }
    {
        Input_Item item("ndy");
        item.annotation = "number of points along y axis for FFT smooth grid";
        read_bcast_int(ndy);
        this->add_item(item);
    }
    {
        Input_Item item("ndz");
        item.annotation = "number of points along z axis for FFT smooth grid";
        read_bcast_int(ndz);
        this->add_item(item);
    }
    {
        Input_Item item("cell_factor");
        item.annotation = "used in the construction of the pseudopotential tables";
        read_bcast_double(cell_factor);
        this->add_item(item);
    }
    {
        Input_Item item("pw_seed");
        item.annotation = "random seed for initializing wave functions";
        read_bcast_int(pw_seed);
        this->add_item(item);
    }
}
} // namespace ModuleIO
