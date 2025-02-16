// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/10/05
// Ref: PHYSICAL REVIEW RESEARCH 5, 013036 (2023)
// This file tests the NCLibxc library.
// how to compile as an independent program:
/*
source  /public1/home/scg0213/software-scg0213/libxc-7.0.0/install/libxc.sh

g++ -std=c++17 -o my_program NCLibxc.cpp LebedevGrid.cpp interface_to_libxc.cpp test_NCLibxc.cpp FibonacciGrid.cpp -I/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/include -L/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/lib -lxc
*/

#include "NCLibxc.h"
#include <iostream>
#include <vector>
#include <array>
#include <complex>

int main() {
        NCLibxc::print_NCLibxc();
        NCLibxc::lda_mc_lc_test();
        NCLibxc::lda_mc_collinear_test();
        NCLibxc::gga_collinear_test();
}