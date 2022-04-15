#ifndef TOOL_CHECK_H
#define TOOL_CHECK_H

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
/**
 * @brief Check the next input from ifs is std::string.
 * This is a global function.
 *
 * @param[in] ifs The input file stream
 * @param[in] name_in The name for checking
 * @param[in] quit Whether call WARNING_QUIT to quit or not
 */
void CHECK_NAME(std::ifstream &ifs, const std::string &name_in, bool quit = true);

/**
 * @brief Check the next input from ifs is integer.
 * This is a global function.
 *
 * @param[in] ifs The input file stream
 * @param[in] v The int variable for checking
 * @param[in] quit Whether call WARNING_QUIT to quit or not
 */
void CHECK_INT(std::ifstream &ifs, const int &v, bool quit = true);

/**
 * @brief Check the next input from ifs is double.
 * This is a global function.
 *
 * @param[in] ifs The input file stream
 * @param[in] v The double variable for checking
 * @param[in] quit Whether call WARNING_QUIT to quit or not
 */
void CHECK_DOUBLE(std::ifstream &ifs, const double &v, bool quit = true);

/**
 * @brief Check the next input from ifs is std::string.
 * This is a global function.
 *
 * @param[in] ifs The input file stream
 * @param[in] v The std::string variable for checking
 * @param[in] quit Whether call WARNING_QUIT to quit or not
 */
void CHECK_STRING(std::ifstream &ifs, const std::string &v, bool quit = true);
} // namespace ModuleBase

#endif
