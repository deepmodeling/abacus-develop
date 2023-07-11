#ifndef __TDDFTTEST
#define __TDDFTTEST
#include <iostream>

#include "gtest/gtest.h"
#include "module_basis/module_ao/parallel_orbitals.h"

extern int myprow, nprow, ictxt, mypcol, npcol;

class TDDFTTEST : public testing::Test
{
  public:
    static void SetUpTestCase()
    {
        std::cout << "\033[32m"
             << "============================"
             << "\033[0m" << std::endl;
        std::cout << "\033[32m"
             << "=    TDDFT MODULE TEST START  ="
             << "\033[0m" << std::endl;
        std::cout << "\033[32m"
             << "============================"
             << "\033[0m" << std::endl;
    }
    static void TearDownTestCase()
    {
        std::cout << "\033[32m"
             << "============================"
             << "\033[0m" << std::endl;
        std::cout << "\033[32m"
             << "=     TDDFT MODULE TEST END   ="
             << "\033[0m" << std::endl;
        std::cout << "\033[32m"
             << "============================"
             << "\033[0m" << std::endl;
    }
    void SetUp()
    {
        std::cout << "\033[32m"
             << "[    CASE  ]"
             << "\033[0m"
             << " ";
    }
    void TearDown()
    {
    }
};

#endif
