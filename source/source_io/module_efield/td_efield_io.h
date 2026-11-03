#ifndef TD_EFIELD_IO_H
#define TD_EFIELD_IO_H

#include <cstddef>
#include <string>

namespace elecstate
{
class TDFieldManager;
}

namespace ModuleIO
{

/**
 * @brief Prepare electric-field output files for an RT-TDDFT calculation.
 *
 * @param output_dir Output directory including its trailing path separator.
 * @param field_count Number of configured time-dependent electric fields.
 * @param restart Whether the calculation continues from an MD/RT-TDDFT restart.
 */
void prep_td_field_out(const std::string& output_dir, std::size_t field_count, bool restart);

/** @brief Prepare total endpoint E and, for analytic fields, the configured waveform files. */
void prep_td_field_out(const std::string& output_dir, const elecstate::TDFieldManager& manager, bool restart);

/**
 * @brief Append the current electric-field samples to their output files.
 * Times are written in fs with six decimal places, starting at physical time zero.
 * Field components use scientific notation with twelve decimal places; columns are right-aligned.
 *
 * @param manager Time-dependent electric-field state for the current electronic step.
 * @param output_dir Output directory including its trailing path separator.
 */
void write_td_field(const elecstate::TDFieldManager& manager, const std::string& output_dir);

} // namespace ModuleIO

#endif
