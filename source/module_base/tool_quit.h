#ifndef QUIT_H
#define QUIT_H

#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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
 * @brief Close .log files and exit
 *
 */
void QUIT(void);
void QUIT(int ret);

/**
 * @brief Combine the functions of WARNING and QUIT
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
void WARNING_QUIT(const std::string &file, const std::string &description);

/**
 * @brief Check, if true, WARNING_QUIT
 *
 * @param file The file where warning happens
 * @param description The warning information
 */
void CHECK_WARNING_QUIT(const bool error, const std::string &file,const std::string &description);

} // namespace ModuleBase

#endif
