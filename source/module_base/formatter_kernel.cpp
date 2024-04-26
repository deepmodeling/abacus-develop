#include "formatter_kernel.h"
#include <string>
#include <iostream>
#include <vector>


int main()
{
    typedef ABACUSTable t;

    t table(3, 3);
    std::vector<int> col1 = {1, 2, 3};
    std::vector<std::string> col2 = {"a", "b", "c"};
    
    table.fix_fmt("%10d", "%10s");
    table<<"title1"<<"title2";
    table<<col1<<col2;
    std::cout << table.str() << std::endl;
    return 0;
}