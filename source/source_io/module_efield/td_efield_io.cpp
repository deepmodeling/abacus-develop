#include "td_efield_io.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_estate/module_pot/td_field_manager.h"

#include <fstream>
#include <sstream>

namespace
{

std::string td_field_output_path(const std::string& output_dir, const std::size_t field_index)
{
    std::stringstream path;
    path << output_dir << "efield_" << field_index + 1 << ".txt";
    return path.str();
}

} // namespace

namespace ModuleIO
{

void prepare_td_field_output(const std::string& output_dir, const std::size_t field_count, const bool restart)
{
    // A restart continues an existing time series. A fresh calculation opens
    // every configured file with `out`, which truncates any previous series.
    if (restart)
    {
        return;
    }

    for (std::size_t field_index = 0; field_index < field_count; ++field_index)
    {
        const std::string output_path = td_field_output_path(output_dir, field_index);
        std::ofstream output(output_path.c_str(), std::ofstream::out);
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::prepare_td_field_output", "Cannot prepare electric-field file " + output_path + "!");
        }
    }
}

void write_td_field_values(const elecstate::TDFieldManager& manager, const std::string& output_dir)
{
    if (!manager.active())
    {
        return;
    }

    const std::vector<double>& field_values = manager.field_values();
    for (std::size_t field_index = 0; field_index < field_values.size(); ++field_index)
    {
        // Keep one file per input occurrence even when directions repeat.
        const std::string output_path = td_field_output_path(output_dir, field_index);
        std::ofstream output(output_path.c_str(), std::ofstream::app);
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_td_field_values", "Cannot append electric-field file " + output_path + "!");
        }
        // Convert only at the user-visible output boundary: time to fs and the
        // electric field to V/Angstrom.
        output << manager.current_step() * manager.dt() * ModuleBase::AU_to_FS << "\t"
               << field_values[field_index] * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << std::endl;
    }
}

} // namespace ModuleIO
