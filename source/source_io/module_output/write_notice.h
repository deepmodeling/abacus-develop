#ifndef MODULEIO_WRITE_NOTICE_H
#define MODULEIO_WRITE_NOTICE_H

#include <fstream>
#include <string>

#include "spin_tag.h"

namespace ModuleIO
{

/// @brief Announce a written output file in the running log.
///
/// All file-written notices share this single template so that the verb,
/// the preposition and the spin label stay identical across writers:
/// " Write <property><spin> to file: <file_name>".
///
/// @param ofs_running running log stream, open on the writing rank only.
/// @param property name of the written quantity, e.g. "Charge density".
/// @param file_name path of the written file.
inline void write_notice(std::ofstream& ofs_running,
                         const std::string& property,
                         const std::string& file_name)
{
    ofs_running << " Write " << property << " to file: " << file_name << std::endl;
}

/// @brief Announce a spin-resolved output file in the running log.
///
/// @param ofs_running running log stream, open on the writing rank only.
/// @param property name of the written quantity, e.g. "DOS".
/// @param is spin channel index.
/// @param nspin number of spin channels; the label is empty unless nspin == 2.
/// @param file_name path of the written file.
inline void write_notice(std::ofstream& ofs_running,
                         const std::string& property,
                         const int is,
                         const int nspin,
                         const std::string& file_name)
{
    ofs_running << " Write " << property << make_spin_tag(is, nspin) << " to file: " << file_name << std::endl;
}

/// @brief Announce a written NAO-basis matrix file in the running log.
///
/// @param ofs_running running log stream, open on the writing rank only.
/// @param property name of the matrix, e.g. "DM(k)" or "dS/dR (x component)".
/// @param file_name path of the written file.
inline void write_matrix_notice(std::ofstream& ofs_running,
                                const std::string& property,
                                const std::string& file_name)
{
    ofs_running << " Write " << property << " matrix in NAO basis to file: " << file_name << std::endl;
}

/// @brief Announce a written spin-resolved NAO-basis matrix file in the log.
///
/// @param ofs_running running log stream, open on the writing rank only.
/// @param property name of the matrix, e.g. "DM(R)".
/// @param is spin channel index.
/// @param nspin number of spin channels; the label is empty unless nspin == 2.
/// @param file_name path of the written file.
inline void write_matrix_notice(std::ofstream& ofs_running,
                                const std::string& property,
                                const int is,
                                const int nspin,
                                const std::string& file_name)
{
    ofs_running << " Write " << property << make_spin_tag(is, nspin)
                << " matrix in NAO basis to file: " << file_name << std::endl;
}

} // namespace ModuleIO

#endif
