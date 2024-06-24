#include "read_input.h"
#include "read_input_tool.h"


namespace ModuleIO
{
void ReadInput::item_lcao()
{
    {
        Input_Item item("basis_type");
        item.annotation = "PW; LCAO in pw; LCAO";
        read_sync_string(basis_type);
        this->add_item(item);
    }
    {
        Input_Item item("gamma_only");
        item.annotation = "Only for localized orbitals set and gamma point. If set to 1, a fast algorithm is used";
        read_sync_bool(gamma_only);
        this->add_item(item);
    }
    {
        Input_Item item("search_radius");
        item.annotation = "input search radius (Bohr)";
        read_sync_double(search_radius);
        this->add_item(item);
    }
    {
        Input_Item item("search_pbc");
        item.annotation = "input periodic boundary condition";
        read_sync_bool(search_pbc);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_ecut");
        item.annotation = "energy cutoff for LCAO";
        read_sync_double(lcao_ecut);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dk");
        item.annotation = "delta k for 1D integration in LCAO";
        read_sync_double(lcao_dk);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_dr");
        item.annotation = "delta r for 1D integration in LCAO";
        read_sync_double(lcao_dr);
        this->add_item(item);
    }
    {
        Input_Item item("lcao_rmax");
        item.annotation = "max R for 1D two-center integration table";
        read_sync_double(lcao_rmax);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs");
        item.annotation = "output H and S matrix (with precision 8)";
        item.readvalue = [](const Input_Item& item, Parameter& para) {
            int count = item.str_values.size();
            if (count == 1)
            {
                para.input.out_mat_hs[0] = convertstr<int>(item.str_values[0]);
                para.input.out_mat_hs[1] = 8;
            }
            else if (count == 2)
            {
                para.input.out_mat_hs[0] = convertstr<int>(item.str_values[0]);
                para.input.out_mat_hs[1] = convertstr<int>(item.str_values[1]);
            }
            else
            {
                throw std::runtime_error("out_mat_hs should have 1 or 2 values");
            }
        };
        sync_intvec(out_mat_hs, 2);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_hs2");
        item.annotation = "output H(R) and S(R) matrix";
        read_sync_bool(out_mat_hs2);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_dh");
        item.annotation = "output of derivative of H(R) matrix";
        read_sync_bool(out_mat_dh);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_xc");
        item.annotation = "output exchange-correlation matrix in KS-orbital representation";
        read_sync_bool(out_mat_xc);
        this->add_item(item);
    }
    {
        Input_Item item("out_hr_npz");
        item.annotation = "output hr(I0,JR) submatrices in npz format";
        read_sync_bool(out_hr_npz);
        this->add_item(item);
    }
    {
        Input_Item item("out_dm_npz");
        item.annotation = "output dmr(I0,JR) submatrices in npz format";
        read_sync_bool(out_dm_npz);
        this->add_item(item);
    }
    {
        Input_Item item("dm_to_rho");
        item.annotation = "reads dmr in npz format and calculates electron density";
        read_sync_bool(dm_to_rho);
        this->add_item(item);
    }
    {
        Input_Item item("out_interval");
        item.annotation = "interval for printing H(R) and S(R) matrix during MD";
        read_sync_int(out_interval);
        this->add_item(item);
    }
    {
        Input_Item item("out_app_flag");
        item.annotation = "whether output r(R), H(R), S(R), T(R), and dH(R) matrices in an append manner during MD";
        read_sync_bool(out_app_flag);
        this->add_item(item);
    }
    {
        Input_Item item("out_ndigits");
        item.annotation = "the length of decimal part of output data";
        read_sync_int(out_ndigits);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_t");
        item.annotation = "output T(R) matrix";
        read_sync_bool(out_mat_t);
        this->add_item(item);
    }
    {
        Input_Item item("out_element_info");
        item.annotation = "output (projected) wavefunction of each element";
        read_sync_bool(out_element_info);
        this->add_item(item);
    }
    {
        Input_Item item("out_mat_r");
        item.annotation = "output r(R) matrix";
        read_sync_bool(out_mat_r);
        this->add_item(item);
    }
    {
        Input_Item item("out_wfc_lcao");
        item.annotation = "ouput LCAO wave functions, 0, no output 1: text, 2: binary";
        read_sync_int(out_wfc_lcao);
        this->add_item(item);
    }
    {
        Input_Item item("bx");
        item.annotation = "division of an element grid in FFT grid along x";
        read_sync_int(bx);
        this->add_item(item);
    }
    {
        Input_Item item("by");
        item.annotation = "division of an element grid in FFT grid along y";
        read_sync_int(by);
        this->add_item(item);
    }
    {
        Input_Item item("bz");
        item.annotation = "division of an element grid in FFT grid along z";
        read_sync_int(bz);
        this->add_item(item);
    }
    {
        Input_Item item("num_stream");
        item.annotation = "the nstream in compute the LCAO with CUDA";
        read_sync_int(nstream);
        this->add_item(item);
    }
}
} // namespace ModuleIO