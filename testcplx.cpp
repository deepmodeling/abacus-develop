#include <iostream>
#include <complex>
#include "./source/module_base/formatter.h"
int main()
{
    std::complex<double> a = {1.0, 2.0};
    std::cout << FmtCore::format("(%20.10f, %20.10f)", a) << std::endl;
    std::complex<double> b = {3.0, 4.0};
    std::cout << FmtCore::format("(%20.10f, %20.10f) and (%20.10f, %20.10f)", a, b) << std::endl;
    return 0;    
}