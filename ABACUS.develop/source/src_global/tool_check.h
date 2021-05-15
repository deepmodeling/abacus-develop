#ifndef TOOL_CHECK_H 
#define TOOL_CHECK_H 

#include <vector>
#include <valarray>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <cassert>

using namespace std;

//==========================================================
// GLOBAL FUNCTIONS :
// Check the input variables
//==========================================================
void CHECK_NAME(ifstream &ifs, const string &name_in, bool quit=true);
void CHECK_INT(ifstream &ifs, const int &v, bool quit=true);
void CHECK_DOUBLE(ifstream &ifs, const double &v, bool quit=true);
void CHECK_STRING(ifstream &ifs, const string &v, bool quit=true);

#endif
