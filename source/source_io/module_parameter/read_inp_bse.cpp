#include "source_base/tool_quit.h"
#include "read_input.h"
#include "read_input_tool.h"

#include <algorithm>
#include <cctype>

namespace ModuleIO
{
void ReadInput::item_bse()
{
    {
        Input_Item item("bse_tda");
        item.annotation = "whether Tamm-Dancoff Approximation is used (can be 'tda', 'full' or 'both')";
        item.category = "Bethe-Salpeter Equation";
        item.type = "String";
        item.description = "Whether the Tamm-Dancoff approximation is used: 'tda', 'full', or 'both'.";
        item.default_value = "tda";
        read_sync_string(input.bse_tda);
        this->add_item(item);
    }
    {
        Input_Item item("bse_spin_types");
        item.annotation = "which spin type is calculated (can be 'singlet', 'triplet', also for test 'rpa', 'ipa')";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Vector of String (>=1 values)";
        item.description
            = "Spin types for a closed-shell calculation in one task: 'singlet', 'triplet', and the test modes "
              "'rpa' and 'ipa'.";
        item.default_value = "singlet triplet";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            size_t count = item.get_size();
            auto& ist = para.input.bse_spin_types;
            ist.clear();
            for (int i = 0; i < count; i++) { ist.push_back(item.str_values[i]); }
            };
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.bse_spin_types.empty()) { para.input.bse_spin_types={"singlet","triplet"}; }
            };
        sync_stringvec(input.bse_spin_types, para.input.bse_spin_types.size(), "singlet");
        this->add_item(item);
    }
    {
        Input_Item item("bse_mem_save");
        item.annotation = "whether to save memory by adding V and W to BSE matrix directly";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Boolean";
        item.description
            = "Whether to save memory by adding V and W directly to the BSE matrix. When enabled, "
              "bse_ri_hartree is enabled and bse_continue is reset to 0.";
        item.default_value = "false";
        item.reset_value = [](const Input_Item& item, Parameter& para) {
            if (para.input.bse_mem_save == true) { para.input.bse_continue=0; para.input.bse_ri_hartree=true; }
            };
        read_sync_bool(input.bse_mem_save);
        this->add_item(item);
    }
    {
        Input_Item item("bse_ri_hartree");
        item.annotation = "whether to use RI approximation for Hartree term in BSE";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Boolean";
        item.description = "Whether to use the RI approximation for the Hartree term in BSE.";
        item.default_value = "true";
        read_sync_bool(input.bse_ri_hartree);
        this->add_item(item);
    }
    {
        Input_Item item("bse_use_fine_kgrid");
        item.annotation = "Fine k-grid mode for BSE";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Integer";
        item.description
            = "Fine k-grid mode for BSE: 0 uses the coarse k-grid, 1 uses a uniform fine k-grid, and 2 uses a "
              "non-uniform fine k-grid. Modes 1 and 2 require band_kpath_info, band_KS_eigenvector_k_{index}.txt, "
              "KS_band_spin_{index}.txt, and GW_band_spin_{index}.txt.";
        item.default_value = "0";
        read_sync_int(input.bse_use_fine_kgrid);
        this->add_item(item);
    }
    {
        Input_Item item("bse_q_approx_mode");
        item.annotation = "q->kpair mapping mode: 0=exact, 1=coarse q grid, 2=mixed";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Integer";
        item.description
            = "q-to-k-pair mapping mode: 0 uses exact mapping, 1 uses the coarse q-grid approximation, and 2 uses "
              "exact for Γ-close q-points and coarse for other q-points.";
        item.default_value = "0";
        read_sync_int(input.bse_q_approx_mode);
        this->add_item(item);
    }
    {
        Input_Item item("bse_q_approx_threshold");
        item.annotation = "threshold radius (Bohr^-1) for exact q mapping in mode 2";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Real";
        item.description
            = "Threshold radius in Bohr^-1 for exact q-to-k-pair mapping when bse_q_approx_mode is 2.";
        item.default_value = "0.1";
        read_sync_double(input.bse_q_approx_threshold);
        this->add_item(item);
    }
    {
        Input_Item item("out_bse_ab");
        item.annotation = "whether to output the AB matrix to file";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Boolean";
        item.description = "Whether to output the AB matrix to a file.";
        item.default_value = "false";
        read_sync_bool(input.out_bse_ab);
        this->add_item(item);
    }
    {
        Input_Item item("bse_continue");
        item.annotation = "which step to continue from previous BSE calculation";
        item.category = "Bethe-Salpeter Equation";
        item.type = "Integer";
        item.description
            = "Step from which to continue a previous BSE calculation: 0 starts a new calculation; 1 reads A_V; "
              "2 reads A_V and A_W; 3 reads A_V, A_W, and B_V; 4 reads A_V, A_W, B_V, and B_W.";
        item.default_value = "0";
        read_sync_int(input.bse_continue);
        this->add_item(item);
    }
    {
        Input_Item item("plot_istate");
        item.annotation = "which state of exciton to be ploted";
        item.category = "Linear Response TDDFT";
        item.type = "Integer";
        item.description = "The index of the excited state to plot, starting from 0.";
        item.default_value = "0";
        read_sync_int(input.plot_istate);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_plot_type");
        item.annotation = "exciton density plot type: average or conditional";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.description
            = "Exciton density represented when lr_solver is 'plot': 'average' integrates out the other particle, "
              "while 'conditional' fixes one particle at exciton_fixed_coordinate and plots a slice of the other "
              "particle's density.";
        item.default_value = "average";
        read_sync_string(input.exciton_plot_type);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_plot_format");
        item.annotation = "exciton plot format: average supports cube, slice, or both; conditional supports slice only";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.description
            = "The exciton-density output format. Average density supports cube, slice, and both; conditional density "
              "supports slice only.";
        item.default_value = "cube";
        item.unit = "";
        item.set_availability("lr_solver==plot");
        read_sync_string(input.exciton_plot_format);
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const auto lower = [](std::string value) {
                std::transform(value.begin(),
                               value.end(),
                               value.begin(),
                               [](const unsigned char c) { return static_cast<char>(std::tolower(c)); });
                return value;
            };
            if (lower(para.input.lr_solver) != "plot")
            {
                return;
            }
            const std::string format = lower(para.input.exciton_plot_format);
            if (format != "cube" && format != "slice" && format != "both")
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exciton_plot_format must be cube, slice, or both");
            }
            if (lower(para.input.exciton_plot_type) == "conditional" && format != "slice")
            {
                ModuleBase::WARNING_QUIT(
                    "ReadInput", "conditional exciton density only supports exciton_plot_format = slice");
            }
        };
        this->add_item(item);
    }
    {
        Input_Item item("exciton_fixed_coordinate");
        item.annotation = "fixed hole and electron Cartesian coordinates (Bohr): hx hy hz ex ey ez";
        item.category = "Linear Response TDDFT";
        item.type = "Vector of Real (6 values)";
        item.description
            = "Cartesian coordinates in Bohr used by conditional exciton plotting, in the order hole_x hole_y hole_z "
              "electron_x electron_y electron_z.";
        item.default_value = "0.0 0.0 0.0 0.0 0.0 0.0";
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.exciton_fixed_coordinate);
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            if (para.input.exciton_fixed_coordinate.size() != 6)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exciton_fixed_coordinate should have exactly 6 values");
            }
        };
        sync_doublevec(input.exciton_fixed_coordinate, 6, 0.0);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_slice_plane");
        item.annotation = "cross-section plane: ab, bc, or ca";
        item.category = "Linear Response TDDFT";
        item.type = "String";
        item.description
            = "Pair of lattice-vector directions spanning the cross section: 'ab', 'bc', or 'ca'.";
        item.default_value = "ab";
        read_sync_string(input.exciton_slice_plane);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_slice_pos");
        item.annotation = "offset along perpendicular direction (Bohr) for slice";
        item.category = "Linear Response TDDFT";
        item.type = "Real";
        item.description
            = "Offset in Bohr along the remaining lattice-vector direction: c for an ab slice, a for bc, and b for "
              "ca.";
        item.default_value = "0.0";
        read_sync_double(input.exciton_slice_pos);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_slice_npoints");
        item.annotation = "grid points per dimension for slice";
        item.category = "Linear Response TDDFT";
        item.type = "Integer";
        item.description
            = "Target in-plane grid resolution. The final grid uses a uniform number of points per primitive cell "
              "over `exciton_slice_range` and includes both endpoints.";
        item.default_value = "200";
        read_sync_int(input.exciton_slice_npoints);
        this->add_item(item);
    }
    {
        Input_Item item("exciton_slice_range");
        item.annotation = "slice cell range: ustart uend vstart vend";
        item.category = "Linear Response TDDFT";
        item.type = "Vector of Integer (4 values)";
        item.description
            = "The in-plane primitive-cell range of an exciton slice: ustart uend vstart vend. The end values are "
              "exclusive cell boundaries, while grid data include both range endpoints.";
        item.default_value = "-1 2 -1 2";
        item.unit = "primitive cells";
        item.set_availability("lr_solver==plot and exciton_plot_format in [slice, both]");
        item.read_value = [](const Input_Item& item, Parameter& para) {
            parse_expression(item.str_values, para.input.exciton_slice_range);
        };
        item.check_value = [](const Input_Item&, const Parameter& para) {
            const std::vector<int>& range = para.input.exciton_slice_range;
            if (range.size() != 4)
            {
                ModuleBase::WARNING_QUIT("ReadInput", "exciton_slice_range should have exactly 4 integer values");
            }
            else if (range[0] >= range[1] || range[2] >= range[3])
            {
                ModuleBase::WARNING_QUIT("ReadInput",
                                         "exciton_slice_range should satisfy ustart < uend and vstart < vend");
            }
        };
        sync_intvec(input.exciton_slice_range, para.input.exciton_slice_range.size(), 0);
        this->add_item(item);
    }
}
} // namespace ModuleIO
