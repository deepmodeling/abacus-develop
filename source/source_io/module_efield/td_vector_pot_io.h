#ifndef TD_VECTOR_POT_IO_H
#define TD_VECTOR_POT_IO_H

#include "source_base/vector3.h"

#include <string>
#include <vector>

namespace ModuleIO
{

/**
 * @brief Read Cartesian vector potentials for an RT-TDDFT calculation.
 *
 * @param input_dir Input directory including its trailing path separator.
 * @return Vector potentials in atomic units, ordered by file row.
 */
std::vector<ModuleBase::Vector3<double>> read_td_vector_pot(const std::string& input_dir);

/**
 * @brief Prepare the vector-potential output file for an RT-TDDFT calculation.
 *
 * @param output_dir Output directory including its trailing path separator.
 * @param restart Whether the calculation continues from an MD/RT-TDDFT restart.
 */
void prepare_td_vector_pot_output(const std::string& output_dir, bool restart);

/**
 * @brief Append one Cartesian vector-potential sample.
 *
 * @param output_dir Output directory including its trailing path separator.
 * @param electronic_step Zero-based electronic-step index.
 * @param vector_pot Cartesian vector potential in atomic units.
 */
void write_td_vector_pot(const std::string& output_dir, int electronic_step, const ModuleBase::Vector3<double>& vector_pot);

} // namespace ModuleIO

#endif
