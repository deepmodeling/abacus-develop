#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_parameter/parameter.h"
#include <cmath>
#include <iostream>

int main()
{
    int checked = 0;
    for (int window : {1, 3, 5})
    {
        for (int restart : {0, 4})
        {
            for (double threshold : {0.0, -0.01, -0.3})
            {
                for (double slope : {-0.5, 0.0, 0.25})
                {
                    Charge_Mixing mixing;
                    mixing.mixing_restart_last = restart;
                    for (int step = 1; step <= 10; ++step)
                    {
                        const double density_error = std::pow(10.0, slope * step);
                        const bool actual = mixing.if_scf_oscillate(step, density_error, window, threshold);
                        const bool expected = threshold < 0 && step >= window + restart && window > 1 && slope > threshold;
                        if (actual != expected) { return 2; }
                        ++checked;
                    }
                }
            }
        }
    }
    std::cout << "oscillation_cases=" << checked << '\n';
}
