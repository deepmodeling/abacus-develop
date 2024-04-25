#include "formatter_kernel.h"
#include <string>
#include <iostream>
#include <vector>


int main()
{
    typedef ABACUSFormatter f; // ABACUSFormatter is a class defined in formatter_kernel.h

    std::string s = f::format("%d", 1);
    std::cout << s << std::endl;
    s = f::format("%f %f", 1.0, 2.0);
    std::cout << s << std::endl;
    std::string fmtstr = "%d %f %s";
    
    f fmt("%d");
    s = fmt.format(1);
    std::cout << s << std::endl;
    fmt.reset("%20.10f %16.6f");
    s = fmt.format(1.0, 2.0);
    std::cout << s << std::endl;   

    // test the spill-over case
    s = f::format("%1d %1d", 100, 20000000);
    std::cout << s << std::endl;
    s = f::format("%2.1f %2.1f", 1.2345678910, 200.345678910);
    std::cout << s << std::endl;
    return 0;
}