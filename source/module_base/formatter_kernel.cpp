#include "formatter_kernel.h"
#include <string>
#include <iostream>
#include <vector>


int main()
{
    typedef ABACUSTable t;

    t table(3, 4);
    std::vector<int> col1 = {1, 2, 3};
    // const char* will not trigger the bug
    // std::string will trigger the bug
    std::vector<std::string> col2 = {"a123", "b234", "c345"};
    
    table.fix_fmt("%10d", "%10s", "%10d", "%10s");
    table<<"title1"<<"title2"<<"title3"<<"title4";
    table<<col1<<col2<<col1<<col2;
    std::cout << table.str() << std::endl;
    return 0;
}