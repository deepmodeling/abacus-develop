#include "source_io/module_efield/td_efield_io.h"

#include "source_base/constants.h"
#include "source_base/tool_quit.h"
#include "source_estate/module_pot/td_field_manager.h"

#include <fstream>
#include <iomanip>
#include <sstream>

namespace
{

std::string field_out_path(const std::string& output_dir, const std::size_t field_index)
{
    std::stringstream path;
    path << output_dir << "efield_" << field_index + 1 << ".txt";
    return path.str();
}

} // namespace

namespace ModuleIO
{

void prep_td_field_out(const std::string& output_dir, const std::size_t field_count, const bool restart)
{
    // A restart continues an existing time series. A fresh calculation opens
    // every configured file with `out`, which truncates any previous series.
    if (restart)
    {
        return;
    }

    for (std::size_t field_index = 0; field_index < field_count; ++field_index)
    {
        const std::string output_path = field_out_path(output_dir, field_index);
        std::ofstream output(output_path.c_str(), std::ofstream::out);
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::prep_td_field_out", "Cannot prepare electric-field file " + output_path + "!");
        }
    }
}

void prep_td_field_out(const std::string& output_dir, const elecstate::TDFieldManager& manager, const bool restart)
{
    if (!manager.file_source())
    {
        prep_td_field_out(output_dir, manager.fields().size(), restart);
    }
    if (!restart)
    {
        std::ofstream output((output_dir + "efield_total.txt").c_str());
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::prep_td_field_out", "Cannot prepare efield_total.txt.");
        }
        output << std::left << std::setw(10) << "# time(fs)" << std::right;
        for (const char* label: {"Ex(V/Angstrom)", "Ey(V/Angstrom)", "Ez(V/Angstrom)"})
        {
            output << "  " << std::setw(19) << label;
        }
        output << "  # ";
        if (manager.file_source())
        {
            output << "endpoint finite differences\n";
        }
        else
        {
            output << "sum of analytic fields\n";
        }
    }
}

void write_td_field(const elecstate::TDFieldManager& manager, const std::string& output_dir)
{
    const double time_fs = manager.current_step() * manager.dt_ha() * ModuleBase::AU_to_FS;
    // Keep the complete time series, including zero fields outside the analytic window.
    {
        std::ofstream output((output_dir + "efield_total.txt").c_str(), std::ofstream::app);
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_td_field", "Cannot append efield_total.txt.");
        }
        output << std::right << std::fixed << std::setprecision(6) << std::setw(10) << time_fs;
        output << std::scientific << std::setprecision(12);
        for (int direction = 0; direction < 3; ++direction)
        {
            output << "  " << std::setw(19) << manager.efield_ha()[direction] * ModuleBase::Hartree_to_eV / ModuleBase::BOHR_TO_A;
        }
        output << '\n';
    }
    if (manager.file_source())
    {
        return;
    }

    const std::vector<double>& field_values = manager.field_vals_ha();
    for (std::size_t field_index = 0; field_index < field_values.size(); ++field_index)
    {
        // Keep one file per input occurrence even when directions repeat.
        const std::string output_path = field_out_path(output_dir, field_index);
        std::ofstream output(output_path.c_str(), std::ofstream::app);
        if (!output)
        {
            ModuleBase::WARNING_QUIT("ModuleIO::write_td_field", "Cannot append electric-field file " + output_path + "!");
        }
        // Convert only at the user-visible output boundary: time to fs and the
        // electric field to V/Angstrom.
        output << std::right << std::fixed << std::setprecision(6) << std::setw(10) << time_fs << "  " << std::scientific
               << std::setprecision(12) << std::setw(19) << field_values[field_index] * ModuleBase::Hartree_to_eV / ModuleBase::BOHR_TO_A
               << std::endl;
    }
}

} // namespace ModuleIO
