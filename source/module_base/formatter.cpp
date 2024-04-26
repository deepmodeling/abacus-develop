#include "formatter.h"
#include <string>
#include <iostream>
#include <vector>


int main()
{
    std::vector<std::string> titles = {"Title of properties", "Title of values", "Title of values column2", "title4"};
    FmtTable table(titles, 3, {"%10d", "%10s", "%20.10f", "%10s"});

    std::vector<int> col1 = {1, 2, 3};
    std::vector<std::string> col2 = {"This is the first line", "b234", "c345"};
    std::vector<double> col3 = {1.1, 2.2, 3.3};
    std::vector<std::string> col4 = {"a123", "This is the second line", "c345"};
    
    table<<col1<<col2<<col3<<col4;

    std::cout << table.str() << std::endl;
    return 0;
}