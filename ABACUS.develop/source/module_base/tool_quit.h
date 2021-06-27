#ifndef QUIT_H 
#define QUIT_H 

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
// GLOBAL FUNCTION :
// NAME : WARNING( write information into ofs_warning)
// NAME : QUIT( exit the running program)
// NAME : WARNING_QUIT( write information into
// 		  ofs_warning , and then quit)
//==========================================================
void WARNING(const string &file,const string &description);
void QUIT(void);
void WARNING_QUIT(const string &file,const string &description);


#endif
