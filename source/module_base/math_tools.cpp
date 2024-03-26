#include "module_base/math_tools.h"
#include <cassert>

namespace ModuleBase
{

std::pair<double, double> linear_fit(const std::vector<double>& x, const std::vector<double>& y)
{
    assert(x.size() == y.size()); // the size of x and y should be the same
    assert(x.size() > 0); // the size of x should be larger than 0
    // get the size
    const int n = x.size();

  // calculate the sum of x, y, xy, xx
    double x_sum = 0.0, y_sum = 0.0, x2_sum = 0.0, xy_sum = 0.0;
    for(int i = 0; i < n; i++) 
    {
        x_sum += x[i];
        y_sum += y[i];
        x2_sum += x[i] * x[i];
        xy_sum += x[i] * y[i];
    }

    double x_mean = x_sum / n;
    double y_mean = y_sum / n;

    //
    double slope = (n * xy_sum - x_sum * y_sum) / (n * x2_sum - x_sum * x_sum);
    double b = y_mean - slope * x_mean;

    return {slope, b};
}

} // namespace ModuleBase
