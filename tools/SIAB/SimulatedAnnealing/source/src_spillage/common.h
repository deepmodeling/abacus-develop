#ifndef COMMON_H
#define COMMON_H

#include <cstdlib>
#include <new>
#include <cassert>

#include <complex>
#include <cmath>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <vector>

#include "../src_tools/mathzone.h"
#include "../src_tools/realarray.h"
#include "../src_tools/intarray.h"
#include "../src_tools/matrix.h"
#include "../src_tools/matrix3.h"
#include "../src_tools/complexmatrix.h"
#include "../src_tools/complexarray.h"
#include "../src_tools/lapack_connector.h"
#include "../src_tools/Random.h"
#include "../src_tools/timer.h" // added 2009-4-17
#include "../src_tools/inverse_matrix.h"

#ifdef __MPI
#include <mpi.h>
#endif

template <class T>
static void READ_VALUE(std::ifstream& ifs, T& v) {
    ifs >> v;
    ifs.ignore(150, '\n');
    return;
}

bool SCAN_BEGIN(std::ifstream& ifs, const std::string& TargetName, const bool restart = 1, const bool quit = true);
void SCAN_END(std::ifstream& ifs, const std::string& TargetName);
void TITLE(const std::string& class_name, const std::string& function_name);
void TITLE(std::ofstream& ofs, const std::string& class_name, const std::string& function_name);
void PRINTCM(const std::string& s, const ComplexMatrix& m);
//==========================================================
// GLOBAL FUNCTION :
// NAME : ZEROS
// set elements of u as zero which u is 1_d complex array
//==========================================================
template <class T>
void ZEROS(std::complex<T>* u, const int n) {
    assert(n != 0);
    assert(u != 0);
    for (int i = 0; i < n; i++) {
        u[i] = std::complex<T>(0.0, 0.0);
    }
    return;
}

template <class T>
void ZEROS(T* u, const int n) {
    assert(n >= 0);
    for (int i = 0; i < n; i++) {
        u[i] = 0;
    }
}
void WARNING_QUIT(const std::string& file, const std::string& description);
void BLOCK_HERE(const std::string& description);
void QUIT(); // mohan add 2010-06-12

template <class T>
void OUT(std::ofstream& ofs, const std::string& name, const T& a) {
    std::stringstream name2;
    name2 << "[" << name << "]";
    ofs << " " << std::setw(40) << name2.str() << " = " << a << std::endl;
    //    ofs << " " << name << a << std::endl;
    return;
}

template <class T>
void OUT(std::ofstream& ofs, const std::string& name, const T& x, const T& y, const T& z) {
    std::stringstream name2;
    name2 << "[" << name << "]";
    ofs << " " << std::setw(40) << name2.str() << " = " << x << ", " << y << ", " << z << std::endl;
    return;
}

#endif
