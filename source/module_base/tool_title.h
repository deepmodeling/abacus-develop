#ifndef TITLE_H 
#define TITLE_H 

#include <vector>
#include <valarray>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <cassert>

namespace ModuleBase
{

void TITLE(const std::string &class_function_name,const bool disable=true);
void TITLE(const std::string &class_name,const std::string &function_name,const bool disable=true);
void TITLE(std::ofstream &ofs,const std::string &class_name,const std::string &function_name,const bool disable=true);
}

#endif
