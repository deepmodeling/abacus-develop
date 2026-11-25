#include "source_io/module_efield/td_vecpot_io.h"

#include "source_base/tool_quit.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

namespace
{

const char* VECTOR_POT_FILENAME = "vector_pot.txt";

std::string vecpot_path(const std::string& directory)
{
    return directory + VECTOR_POT_FILENAME;
}

void wr_vecpot_header(std::ofstream& output)
{
    output << std::left << std::setw(6) << "#istep" << std::right;
    for (const char* label: {"A_x", "A_y", "A_z"})
    {
        output << "  " << std::setw(19) << label;
    }
    output << "  # endpoint values, Hartree atomic units" << std::endl;
}

bool is_blank_or_comment(const std::string& line)
{
    const std::string::size_type first = line.find_first_not_of(" \t\r\n");
    return first == std::string::npos || line[first] == '#';
}

} // namespace

namespace ModuleIO
{

std::vector<ModuleBase::Vector3<double>> read_td_vecpot(const std::string& input_dir)
{
    const std::string input_path = vecpot_path(input_dir);
    std::ifstream input(input_path.c_str());
    if (!input)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_td_vecpot", "Cannot open vector-potential file " + input_path + "!");
    }

    std::vector<ModuleBase::Vector3<double>> vector_potentials;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        if (is_blank_or_comment(line))
        {
            continue;
        }

        std::istringstream row(line);
        int step_label = 0;
        ModuleBase::Vector3<double> vector_pot;
        if (!(row >> step_label >> vector_pot[0] >> vector_pot[1] >> vector_pot[2]))
        {
            ModuleBase::WARNING_QUIT("ModuleIO::read_td_vecpot",
                                     "Invalid vector-potential data on line " + std::to_string(line_number) + " of " + input_path + "!");
        }
        row >> std::ws;
        for (int direction = 0; direction < 3; ++direction)
        {
            if (!std::isfinite(vector_pot[direction]))
            {
                ModuleBase::WARNING_QUIT("ModuleIO::read_td_vecpot", "Non-finite vector potential in " + input_path);
            }
        }
        if (!row.eof())
        {
            ModuleBase::WARNING_QUIT("ModuleIO::read_td_vecpot",
                                     "Unexpected vector-potential data on line " + std::to_string(line_number) + " of " + input_path + "!");
        }
        vector_potentials.push_back(vector_pot);
    }

    if (vector_potentials.empty())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_td_vecpot", "No vector-potential data found in " + input_path + "!");
    }
    return vector_potentials;
}

void check_td_vecpot_paths(const std::string& input_dir, const std::string& output_dir)
{
    struct stat input_state;
    struct stat output_state;
    const std::string input_path = vecpot_path(input_dir);
    const std::string output_path = vecpot_path(output_dir);
    // Different paths can name the same file through symbolic or hard links; compare file identity before writing.
    if (stat(input_path.c_str(), &input_state) == 0 && stat(output_path.c_str(), &output_state) == 0
        && input_state.st_dev == output_state.st_dev && input_state.st_ino == output_state.st_ino)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::check_td_vecpot_paths", "Vector-potential output would overwrite its input file.");
    }
}

void prep_td_vecpot_out(const std::string& output_dir, const bool restart)
{
    const std::string output_path = vecpot_path(output_dir);
    if (restart)
    {
        std::ifstream existing(output_path.c_str(), std::ifstream::binary | std::ifstream::ate);
        if (existing && existing.tellg() > 0)
        {
            return;
        }
    }

    std::ofstream output(output_path.c_str(), std::ofstream::out);
    if (!output)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::prep_td_vecpot_out", "Cannot prepare vector-potential file " + output_path + "!");
    }
    wr_vecpot_header(output);
}

void write_td_vecpot(const std::string& output_dir, const int electronic_step, const ModuleBase::Vector3<double>& vector_pot)
{
    const std::string output_path = vecpot_path(output_dir);
    std::ofstream output(output_path.c_str(), std::ofstream::app);
    if (!output)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_td_vecpot", "Cannot append vector-potential file " + output_path + "!");
    }

    output << std::right << std::setw(6) << electronic_step + 1 << std::scientific << std::setprecision(12);
    for (int direction = 0; direction < 3; ++direction)
    {
        output << "  " << std::setw(19) << vector_pot[direction];
    }
    output << std::endl;
}

} // namespace ModuleIO
