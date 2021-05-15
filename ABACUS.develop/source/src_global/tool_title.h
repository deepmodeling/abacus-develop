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

using namespace std;

void TITLE(const string &class_function_name);
void TITLE(const string &class_name,const string &function_name);
void TITLE(ofstream &ofs,const string &class_name,const string &function_name);

#endif
