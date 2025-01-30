//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_CONVERT_H
#define INPUT_CONVERT_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

namespace Input_Conv
{

/**
 * @brief template bridge codes for converting string to other types
 *
 */
void tmp_convert();

/**
 * @brief Pass the data members from the INPUT instance(defined in
 * module_io/input.cpp) to GlobalV and GlobalC.
 */
void Convert();

/**
 * @brief To parse input parameters as expressions into vectors
 *
 * @tparam T
 * @param fn  (string): expressions such as "3*1 0 2*0.5 3*0"
 * @param vec (vector): stores parsing results,
 *            for example, "3*1 0 2*0.5 1*1.5" can be parsed as
 *            [1, 1, 1, 0, 0.5, 0.5, 1.5]
 */
template <typename T>
void parse_expression(const std::string& fn, std::vector<T>& vec)
{
    ModuleBase::TITLE("Input_Conv", "parse_expression");
    int count = 0;

    // Update the regex pattern to handle scientific notation
    std::string pattern("([-+]?[0-9]+\\*[-+]?[0-9.eE+-]+|[-+]?[0-9,.eE+-]+)");

    std::vector<std::string> str;
    std::stringstream ss(fn);
    std::string section;

    // Split the input string into substrings by spaces
    while (ss >> section)
    {
        int index = 0;
        if (str.empty())
        {
            while (index < section.size() && std::isspace(section[index]))
            {
                index++;
            }
        }
        section.erase(0, index);
        str.push_back(section);
    }

    // Compile the regular expression
    regex_t reg;
    regcomp(&reg, pattern.c_str(), REG_EXTENDED);
    regmatch_t pmatch[1];
    const size_t nmatch = 1;

    // Loop over each section and apply regex to extract numbers
    for (size_t i = 0; i < str.size(); ++i)
    {
        if (str[i] == "")
        {
            continue;
        }
        int status = regexec(&reg, str[i].c_str(), nmatch, pmatch, 0);
        std::string sub_str = "";

        // Extract the matched substring
        for (size_t j = pmatch[0].rm_so; j != pmatch[0].rm_eo; ++j)
        {
            sub_str += str[i][j];
        }

        // Check if the substring contains multiplication (e.g., "2*3.14")
        std::string sub_pattern("\\*");
        regex_t sub_reg;
        regcomp(&sub_reg, sub_pattern.c_str(), REG_EXTENDED);
        regmatch_t sub_pmatch[1];
        const size_t sub_nmatch = 1;

        if (regexec(&sub_reg, sub_str.c_str(), sub_nmatch, sub_pmatch, 0) == 0)
        {
            size_t pos = sub_str.find("*");
            int num = stoi(sub_str.substr(0, pos));
            assert(num >= 0);
            T occ = stof(sub_str.substr(pos + 1, sub_str.size()));
            
            // Add the value to the vector `num` times
            for (size_t k = 0; k != num; k++)
            {
                vec.emplace_back(occ);
            }
        }
        else
        {
            // Handle scientific notation and convert to T
            std::stringstream convert;
            convert << sub_str;
            T occ;
            convert >> occ;
            vec.emplace_back(occ);
        }

        regfree(&sub_reg);
    }

    regfree(&reg);
}

#ifdef __LCAO
/**
 * @brief convert units of different parameters
 *
 * @param params input parameter
 * @param c coefficients of unit conversion
 * @return parame*c : parameter after unit vonversion
 */
std::vector<double> convert_units(std::string params, double c);

/**
 * @brief read paramers of electric field for tddft and convert units
 */
void read_td_efield();
#endif

} // namespace Input_Conv

#endif // Input_Convert
