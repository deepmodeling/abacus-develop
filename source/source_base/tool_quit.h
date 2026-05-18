#ifndef QUIT_H
#define QUIT_H

#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

namespace ModuleBase
{

//==========================================================
// GLOBAL FUNCTION :
// NAME : WARNING( write information into GlobalV::ofs_warning)
// NAME : QUIT( exit the running program)
// NAME : WARNING_QUIT( write information into
// 		  GlobalV::ofs_warning , and then quit)
//==========================================================
/**
 * @brief Print out warning information in warning.log file
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
void WARNING(const std::string &file, const std::string &description);

/**
 * @brief Inject the global output directory used by QUIT/WARNING_QUIT/CHECK_WARNING_QUIT
 *        for log path resolution and user-facing messages.
 *
 * Caller-injected (typically once after input parameters are read).
 * If never set, paths fall back to CWD.
 */
void set_quit_out_dir(const std::string& dir);

/**
 * @brief Close .log files and exit
 *
 */
[[noreturn]] void QUIT(void);

/**
 * @brief Close .log files and exit
 *
 */
[[noreturn]] void QUIT(int ret);

/**
 * @brief Combine the functions of WARNING and QUIT
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
[[noreturn]] void WARNING_QUIT(const std::string& file, const std::string& description);

/**
 * @brief Combine the functions of WARNING and QUIT
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
[[noreturn]] void WARNING_QUIT(const std::string& file, const std::string& description, int ret);

/**
 * @brief Check, if true, WARNING_QUIT
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
void CHECK_WARNING_QUIT(const bool error, const std::string &file,const std::string &calculation,const std::string &description);

} // namespace ModuleBase

#endif
