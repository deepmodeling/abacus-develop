#ifndef __TDDFTTEST
#define __TDDFTTEST
#include <iostream>

#include "gtest/gtest.h"
#include "module_basis/module_ao/parallel_orbitals.h"

using namespace std;
extern int myprow, nprow, ictxt, mypcol, npcol;

class TDDFTTEST : public testing::Test
{
  public:
    static void SetUpTestCase()
    {
        cout << "\033[32m"
             << "============================"
             << "\033[0m" << endl;
        cout << "\033[32m"
             << "=    TDDFT MODULE TEST START  ="
             << "\033[0m" << endl;
        cout << "\033[32m"
             << "============================"
             << "\033[0m" << endl;
    }
    static void TearDownTestCase()
    {
        cout << "\033[32m"
             << "============================"
             << "\033[0m" << endl;
        cout << "\033[32m"
             << "=     TDDFT MODULE TEST END   ="
             << "\033[0m" << endl;
        cout << "\033[32m"
             << "============================"
             << "\033[0m" << endl;
    }
    void SetUp()
    {
        cout << "\033[32m"
             << "[    CASE  ]"
             << "\033[0m"
             << " ";
    }
    void TearDown()
    {
    }
};

#endif
