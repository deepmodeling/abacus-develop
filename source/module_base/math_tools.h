#ifndef MATH_TOOLS_H
#define MATH_TOOLS_H

#include <vector>

//========================================================
// some fundamental math tools
//========================================================
namespace ModuleBase
{
    /**
     * @brief linear fit by using least square method
     * 
     * @param x [in] x data
     * @param y [in] y data
     * @return std::pair<double, double> slope and intercept
     */
    std::pair<double, double> linear_fit(const std::vector<double>& x, const std::vector<double>& y);
}

#endif