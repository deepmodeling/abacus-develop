//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_CONVERT_H
#define INPUT_CONVERT_H

#include <regex.h>
#include <stdio.h>
#include <string.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

namespace Input_Conv
{
void Convert(void);

// Function `parse_expression` is used to parse input parameters as expressions into vectors
// fn (string): expressions such as "3*1 0 2*0.5 3*0"
// arr (vector): stores parsing results, for example, "3*1 0 2*0.5 1*1.5" can be parsed as [1, 1, 1, 0, 0.5, 0.5, 1.5]
template <typename T> void parse_expression(const std::string &fn, std::vector<T> &arr);

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

/**
 * @brief struct for storing parameters of non-collinear spin-constrained DFT
 */
struct ScElementData {
    int index;
    std::vector<double> lambda;
    std::vector<double> sc_mag;
    double sc_spin_val;
    double sc_spin_angle1;
    double sc_spin_angle2;
};

/**
 * @brief parse json input file for non-collinear spin-constrained DFT
 */
void parseScJsonFile(const std::string& filename, std::map<std::string, std::vector<ScElementData>>& data);

} // namespace Input_Conv

#endif // Input_Convert
