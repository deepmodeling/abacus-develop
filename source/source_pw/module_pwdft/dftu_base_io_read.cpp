#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace DFTU_BASE
{

const char* const read_occup_context = "DFTU_BASE::read_occup_m";

struct OccMatReadState
{
    OccMatReadState(std::ifstream& input,
                    const std::string& filename,
                    const UnitCell& ucell,
                    OccupationMatrix& occ,
                    int nspin,
                    OccmatSocLayout soc_layout)
        : input(input),
          filename(filename),
          ucell(ucell),
          occ(occ),
          nspin(nspin),
          soc_layout(soc_layout),
          line_number(0)
    {
    }

    std::ifstream& input;
    const std::string& filename;
    const UnitCell& ucell;
    OccupationMatrix& occ;
    int nspin;
    OccmatSocLayout soc_layout;
    int line_number;
};

[[noreturn]] static void read_occup_format_error(const OccMatReadState& state);
static bool read_occup_data_line(OccMatReadState& state, std::string& line);
static void read_occup_validate_configuration(const UnitCell& ucell,
                                               const std::vector<int>& l_channel,
                                               int nspin,
                                               int npol,
                                               OccmatSocLayout soc_layout);
static void read_occup_atom_block(OccMatReadState& state, int type, int iat, int l);

void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const std::string& fn,
                  int nspin,
                  int npol,
                  OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "read_occup_m");
    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    read_occup_validate_configuration(ucell, l_channel, nspin, npol, soc_layout);

    std::ifstream input(fn.c_str());
    if (!input)
    {
        ModuleBase::WARNING_QUIT(read_occup_context,
                                 "Can not open compact occupation-matrix file " + fn);
    }

    GlobalV::ofs_running << " DFT+U: read occupation matrix from " << fn << std::endl;
    std::cout << " DFT+U: read occupation matrix from " << fn << std::endl;

    OccMatReadState state(input, fn, ucell, occ, nspin, soc_layout);
    int atom_count = 0;
    for (int type = 0; type < ucell.ntype; ++type)
    {
        const int l = l_channel[type];
        if (l < 0)
        {
            continue;
        }
        for (int atom = 0; atom < ucell.atoms[type].na; ++atom)
        {
            read_occup_atom_block(state, type, ucell.itia2iat(type, atom), l);
            ++atom_count;
        }
    }

    std::string extra_line;
    if (atom_count == 0 || read_occup_data_line(state, extra_line))
    {
        read_occup_format_error(state);
    }
}

[[noreturn]] static void read_occup_format_error(const OccMatReadState& state)
{
    std::ostringstream message;
    message << "Invalid compact occupation-matrix format in " << state.filename;
    if (state.line_number > 0)
    {
        message << " at line " << state.line_number;
    }
    message << ".";
    ModuleBase::WARNING_QUIT(read_occup_context, message.str());
}

static bool read_occup_data_line(OccMatReadState& state, std::string& line)
{
    while (std::getline(state.input, line))
    {
        ++state.line_number;
        const std::string::size_type first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#')
        {
            continue;
        }
        const std::string::size_type comment = line.find('#', first);
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }
        return true;
    }
    return false;
}

static void read_occup_validate_configuration(
    const UnitCell& ucell,
    const std::vector<int>& l_channel,
    const int nspin,
    const int npol,
    const OccmatSocLayout soc_layout)
{
    bool valid = static_cast<int>(l_channel.size()) == ucell.ntype;
    valid = valid && (nspin == 1 || nspin == 2 || nspin == 4);
    valid = valid && npol == ((nspin == 4) ? 2 : 1);
    valid = valid && (soc_layout == DFTU_BASE::SOC_LAYOUT_PAULI
                      || soc_layout == DFTU_BASE::SOC_LAYOUT_SPIN_BASIS_REAL);
    for (int type = 0; valid && type < ucell.ntype; ++type)
    {
        valid = l_channel[type] >= -1 && l_channel[type] <= ucell.atoms[type].nwl;
    }
    if (!valid)
    {
        ModuleBase::WARNING_QUIT(read_occup_context,
                                 "Invalid DFT+U occupation reader configuration.");
    }
}

static void read_occup_atom_header(OccMatReadState& state,
                                   const int type,
                                   const int iat,
                                   const int l)
{
    std::string line;
    if (!read_occup_data_line(state, line))
    {
        read_occup_format_error(state);
    }

    std::istringstream header(line);
    std::string element;
    std::string atom_keyword;
    std::string l_keyword;
    std::string mag_keyword;
    int atom_index = 0;
    int file_l = -1;
    if (!(header >> element >> atom_keyword >> atom_index >> l_keyword >> file_l
          >> mag_keyword)
        || element != state.ucell.atoms[type].label || atom_keyword != "Atom"
        || atom_index != iat + 1 || l_keyword != "L" || file_l != l
        || mag_keyword != "mag")
    {
        read_occup_format_error(state);
    }

    const int component_count = (state.nspin == 4) ? 3 : 1;
    for (int component = 0; component < component_count; ++component)
    {
        double magnetization = 0.0;
        if (!(header >> magnetization) || !std::isfinite(magnetization))
        {
            read_occup_format_error(state);
        }
    }
    std::string extra;
    if (header >> extra)
    {
        read_occup_format_error(state);
    }
}

static void read_occup_block_label(OccMatReadState& state, const int block)
{
    std::string line;
    if (!read_occup_data_line(state, line))
    {
        read_occup_format_error(state);
    }

    std::istringstream label(line);
    std::string spin_keyword;
    const bool offdiagonal = state.nspin == 4 && (block == 1 || block == 2);
    if (offdiagonal)
    {
        std::string spin_pair;
        std::string component;
        const char* expected_component = (block == 1) ? "re" : "im";
        if (!(label >> spin_keyword >> spin_pair >> component)
            || spin_keyword != "spin" || spin_pair != "12"
            || component != expected_component)
        {
            read_occup_format_error(state);
        }
    }
    else
    {
        int spin = 0;
        std::string nelec_keyword;
        double nelec = 0.0;
        const int expected_spin = (block == 0) ? 1 : 2;
        if (!(label >> spin_keyword >> spin >> nelec_keyword >> nelec)
            || spin_keyword != "spin" || spin != expected_spin
            || nelec_keyword != "nelec" || !std::isfinite(nelec))
        {
            read_occup_format_error(state);
        }
    }

    std::string extra;
    if (label >> extra)
    {
        read_occup_format_error(state);
    }
}

static std::vector<double> read_occup_matrix(OccMatReadState& state,
                                             const int size)
{
    std::vector<double> matrix(size * size, 0.0);
    for (int row = 0; row < size; ++row)
    {
        std::string line;
        if (!read_occup_data_line(state, line))
        {
            read_occup_format_error(state);
        }

        std::istringstream row_stream(line);
        for (int column = 0; column < size; ++column)
        {
            double value = 0.0;
            if (!(row_stream >> value) || !std::isfinite(value))
            {
                read_occup_format_error(state);
            }
            matrix[row * size + column] = value;
        }

        std::string extra;
        if (row_stream >> extra)
        {
            read_occup_format_error(state);
        }
    }
    return matrix;
}

static void read_occup_store_blocks(OccMatReadState& state,
                                    const int iat,
                                    const int l,
                                    const std::vector<std::vector<double>>& blocks)
{
    const int size = 2 * l + 1;
    if (state.nspin != 4)
    {
        for (int spin = 0; spin < 2; ++spin)
        {
            state.occ.set_flat(iat, l, spin, blocks[spin]);
        }
        return;
    }

    const int block_size = size * size;
    std::vector<double> matrix(4 * block_size, 0.0);
    if (state.soc_layout == DFTU_BASE::SOC_LAYOUT_PAULI)
    {
        for (int index = 0; index < block_size; ++index)
        {
            matrix[0 * block_size + index] = blocks[0][index] + blocks[3][index];
            matrix[1 * block_size + index] = 2.0 * blocks[1][index];
            matrix[2 * block_size + index] = 2.0 * blocks[2][index];
            matrix[3 * block_size + index] = blocks[0][index] - blocks[3][index];
        }
    }
    else
    {
        const int matrix_size = 2 * size;
        for (int row = 0; row < size; ++row)
        {
            const int block_offset = row * size;
            std::copy_n(blocks[0].begin() + block_offset,
                        size,
                        matrix.begin() + row * matrix_size);
            std::copy_n(blocks[1].begin() + block_offset,
                        size,
                        matrix.begin() + row * matrix_size + size);
            std::copy_n(blocks[1].begin() + block_offset,
                        size,
                        matrix.begin() + (size + row) * matrix_size);
            std::copy_n(blocks[3].begin() + block_offset,
                        size,
                        matrix.begin() + (size + row) * matrix_size + size);
        }
    }
    state.occ.set_flat(iat, l, 0, matrix);
}

static void read_occup_atom_block(OccMatReadState& state,
                                  const int type,
                                  const int iat,
                                  const int l)
{
    read_occup_atom_header(state, type, iat, l);

    const int size = 2 * l + 1;
    const int block_count = (state.nspin == 4) ? 4 : 2;
    std::vector<std::vector<double>> blocks(block_count);
    for (int block = 0; block < block_count; ++block)
    {
        read_occup_block_label(state, block);
        blocks[block] = read_occup_matrix(state, size);
    }
    read_occup_store_blocks(state, iat, l, blocks);
}

} // namespace DFTU_BASE
