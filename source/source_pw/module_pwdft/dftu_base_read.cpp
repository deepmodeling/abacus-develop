#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/unitcell.h"

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
    OccMatReadState(std::ifstream& input, const std::string& filename)
        : input(input), filename(filename), line_number(0)
    {
    }

    std::ifstream& input;
    const std::string& filename;
    int line_number;
};

[[noreturn]] static void read_occup_fail(const OccMatReadState& state,
                                         const std::string& message)
{
    std::ostringstream os;
    os << "invalid compact occupation-matrix file " << state.filename;
    if (state.line_number > 0)
    {
        os << " at line " << state.line_number;
    }
    os << ": " << message;
    ModuleBase::WARNING_QUIT(read_occup_context, os.str());
}

static bool read_occup_next_data_line(OccMatReadState& state, std::string& line)
{
    while (std::getline(state.input, line))
    {
        ++state.line_number;
        const std::string::size_type first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#')
        {
            continue;
        }
        return true;
    }
    return false;
}

static std::string read_occup_require_data_line(OccMatReadState& state,
                                                const std::string& expected)
{
    std::string line;
    if (!read_occup_next_data_line(state, line))
    {
        read_occup_fail(state, "unexpected end of file while reading " + expected);
    }
    return line;
}

static void read_occup_require_line_end(OccMatReadState& state,
                                        std::istringstream& line_stream,
                                        const std::string& expected)
{
    std::string extra;
    if (line_stream >> extra)
    {
        read_occup_fail(state, "unexpected token '" + extra + "' after " + expected);
    }
}

static double read_occup_finite_double(OccMatReadState& state,
                                       std::istringstream& line_stream,
                                       const std::string& expected)
{
    double value = 0.0;
    if (!(line_stream >> value) || !std::isfinite(value))
    {
        read_occup_fail(state, "expected a finite number for " + expected);
    }
    return value;
}

static void read_occup_atom_header(OccMatReadState& state,
                                   const UnitCell& ucell,
                                   const int expected_type,
                                   const int expected_iat,
                                   const int expected_l,
                                   const int nspin)
{
    const std::string line = read_occup_require_data_line(state, "an atom header");
    std::istringstream line_stream(line);
    std::string element;
    std::string atom_keyword;
    std::string l_keyword;
    std::string mag_keyword;
    int atom_index = 0;
    int l = -1;
    if (!(line_stream >> element >> atom_keyword >> atom_index >> l_keyword >> l
          >> mag_keyword))
    {
        read_occup_fail(state, "expected compact header '<element> Atom <index> L <l> mag ...'");
    }
    if (atom_keyword != "Atom" || l_keyword != "L" || mag_keyword != "mag")
    {
        read_occup_fail(state, "only the compact 'Atom <index> L <l> mag ...' header is supported");
    }
    if (element != ucell.atoms[expected_type].label)
    {
        read_occup_fail(state, "element label '" + element + "' does not match the unit cell");
    }
    if (atom_index != expected_iat + 1)
    {
        read_occup_fail(state, "atom index does not match the unit-cell order");
    }
    if (l != expected_l)
    {
        read_occup_fail(state, "L does not match the configured DFT+U channel");
    }

    const int magnetization_components = (nspin == 4) ? 3 : 1;
    for (int component = 0; component < magnetization_components; ++component)
    {
        read_occup_finite_double(state, line_stream, "magnetization");
    }
    read_occup_require_line_end(state, line_stream, "the compact atom header");
}

static void read_occup_spin_nelec_label(OccMatReadState& state,
                                        const int expected_spin)
{
    const std::string line = read_occup_require_data_line(state, "a spin block label");
    std::istringstream line_stream(line);
    std::string spin_keyword;
    std::string nelec_keyword;
    int spin = 0;
    if (!(line_stream >> spin_keyword >> spin >> nelec_keyword)
        || spin_keyword != "spin" || spin != expected_spin
        || nelec_keyword != "nelec")
    {
        read_occup_fail(state, "expected compact label 'spin "
                                   + std::to_string(expected_spin) + " nelec <value>'");
    }
    read_occup_finite_double(state, line_stream, "nelec");
    read_occup_require_line_end(state, line_stream, "the spin block label");
}

static void read_occup_soc_offdiagonal_label(OccMatReadState& state,
                                              const std::string& component)
{
    const std::string line
        = read_occup_require_data_line(state, "the spin 12 " + component + " label");
    std::istringstream line_stream(line);
    std::string spin_keyword;
    std::string block;
    std::string actual_component;
    if (!(line_stream >> spin_keyword >> block >> actual_component)
        || spin_keyword != "spin" || block != "12"
        || actual_component != component)
    {
        read_occup_fail(state, "expected compact label 'spin 12 " + component + "'");
    }
    read_occup_require_line_end(state, line_stream, "the spin 12 " + component + " label");
}

static std::vector<double> read_occup_matrix(OccMatReadState& state,
                                             const int size)
{
    std::vector<double> matrix(size * size, 0.0);
    for (int row = 0; row < size; ++row)
    {
        const std::string line = read_occup_require_data_line(state, "matrix row");
        std::istringstream line_stream(line);
        for (int column = 0; column < size; ++column)
        {
            const int index = row * size + column;
            matrix[index]
                = read_occup_finite_double(state, line_stream, "matrix element");
        }
        read_occup_require_line_end(state, line_stream, "a matrix row");
    }
    return matrix;
}

static void read_occup_collinear_atom(OccMatReadState& state,
                                      OccupationMatrix& occ,
                                      const int iat,
                                      const int l)
{
    const int size = 2 * l + 1;
    std::vector<std::vector<double>> spin_matrices(2);
    for (int spin = 0; spin < 2; ++spin)
    {
        read_occup_spin_nelec_label(state, spin + 1);
        spin_matrices[spin] = read_occup_matrix(state, size);
    }

    for (int spin = 0; spin < 2; ++spin)
    {
        for (int row = 0; row < size; ++row)
        {
            for (int column = 0; column < size; ++column)
            {
                const int index = row * size + column;
                occ.set(iat, l, spin, row, column, spin_matrices[spin][index]);
            }
        }
    }
}

static void read_occup_store_soc_pauli(OccupationMatrix& occ,
                                       const int iat,
                                       const int l,
                                       const std::vector<std::vector<double>>& blocks)
{
    ModuleBase::matrix& matrix = occ.mat(iat, l, 0);
    const int block_size = static_cast<int>(blocks[0].size());
    for (int index = 0; index < block_size; ++index)
    {
        matrix.c[0 * block_size + index] = blocks[0][index] + blocks[3][index];
        matrix.c[1 * block_size + index] = 2.0 * blocks[1][index];
        matrix.c[2 * block_size + index] = 2.0 * blocks[2][index];
        matrix.c[3 * block_size + index] = blocks[0][index] - blocks[3][index];
    }
}

static void read_occup_store_soc_spin_basis(
    OccupationMatrix& occ,
    const int iat,
    const int l,
    const std::vector<std::vector<double>>& blocks)
{
    const int size = 2 * l + 1;
    for (int row = 0; row < size; ++row)
    {
        for (int column = 0; column < size; ++column)
        {
            const int index = row * size + column;
            occ.set(iat, l, 0, row, column, blocks[0][index]);
            occ.set(iat, l, 0, row, size + column, blocks[1][index]);
            occ.set(iat, l, 0, size + row, column, blocks[1][index]);
            occ.set(iat, l, 0, size + row, size + column, blocks[3][index]);
        }
    }
}

static void read_occup_soc_atom(OccMatReadState& state,
                                OccupationMatrix& occ,
                                const int iat,
                                const int l,
                                const OccmatSocLayout soc_layout)
{
    const int size = 2 * l + 1;
    std::vector<std::vector<double>> blocks(4);
    read_occup_spin_nelec_label(state, 1);
    blocks[0] = read_occup_matrix(state, size);
    read_occup_soc_offdiagonal_label(state, "re");
    blocks[1] = read_occup_matrix(state, size);
    read_occup_soc_offdiagonal_label(state, "im");
    blocks[2] = read_occup_matrix(state, size);
    read_occup_spin_nelec_label(state, 2);
    blocks[3] = read_occup_matrix(state, size);

    if (soc_layout == SOC_LAYOUT_PAULI)
    {
        read_occup_store_soc_pauli(occ, iat, l, blocks);
    }
    else
    {
        read_occup_store_soc_spin_basis(occ, iat, l, blocks);
    }
}

static void read_occup_validate_configuration(
    const OccMatReadState& state,
    const UnitCell& ucell,
    const std::vector<int>& l_channel,
    const int nspin,
    const int npol,
    const OccmatSocLayout soc_layout)
{
    if (static_cast<int>(l_channel.size()) != ucell.ntype)
    {
        read_occup_fail(state, "l_channel size does not match the number of atom types");
    }
    for (int type = 0; type < ucell.ntype; ++type)
    {
        if (l_channel[type] < -1 || l_channel[type] > ucell.atoms[type].nwl)
        {
            read_occup_fail(state, "a configured DFT+U channel is outside the unit-cell range");
        }
    }
    if (nspin != 1 && nspin != 2 && nspin != 4)
    {
        read_occup_fail(state, "nspin must be 1, 2 or 4");
    }
    const int expected_npol = (nspin == 4) ? 2 : 1;
    if (npol != expected_npol)
    {
        read_occup_fail(state, "npol is inconsistent with nspin");
    }
    if (soc_layout != SOC_LAYOUT_PAULI
        && soc_layout != SOC_LAYOUT_SPIN_BASIS_REAL)
    {
        read_occup_fail(state, "unknown SOC occupation-matrix storage layout");
    }
}

static void read_occup_report_open_error(const std::string& filename)
{
    ModuleBase::WARNING_QUIT(read_occup_context,
                             "Can not open compact occupation-matrix file " + filename);
}

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

    std::ifstream input(fn.c_str());
    if (!input)
    {
        read_occup_report_open_error(fn);
    }

    GlobalV::ofs_running << " DFT+U: read occupation matrix from " << fn << std::endl;
    std::cout << " DFT+U: read occupation matrix from " << fn << std::endl;

    OccMatReadState state(input, fn);
    read_occup_validate_configuration(state, ucell, l_channel, nspin, npol, soc_layout);

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
            const int iat = ucell.itia2iat(type, atom);
            read_occup_atom_header(state, ucell, type, iat, l, nspin);
            if (nspin == 4)
            {
                read_occup_soc_atom(state, occ, iat, l, soc_layout);
            }
            else
            {
                read_occup_collinear_atom(state, occ, iat, l);
            }
            ++atom_count;
        }
    }

    if (atom_count == 0)
    {
        read_occup_fail(state, "the unit cell has no configured DFT+U atom blocks");
    }

    std::string extra_line;
    if (read_occup_next_data_line(state, extra_line))
    {
        read_occup_fail(state, "unexpected content after the final atom block");
    }
}

} // namespace DFTU_BASE
