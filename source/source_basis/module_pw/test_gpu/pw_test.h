#ifndef __PWTEST
#define __PWTEST
#include "gtest/gtest.h"
#include <iostream>
#include <string>
extern int nproc_in_pool, rank_in_pool;
extern std::string precision_flag, device_flag;

class PWTEST : public testing::Test {
  public:
    static void SetUpTestCase() {
        if (rank_in_pool == 0) {
            std::cout << "\033[32m" << "============================" << "\033[0m" << std::endl;
            std::cout << "\033[32m" << "=    PW GPU MODULE TEST START  =" << "\033[0m" << std::endl;
            std::cout << "\033[32m" << "============================" << "\033[0m" << std::endl;
        }
    }
    static void TearDownTestCase() {
        if (rank_in_pool == 0) {
            std::cout << "\033[32m" << "============================" << "\033[0m" << std::endl;
            std::cout << "\033[32m" << "=     PW GPU MODULE TEST END   =" << "\033[0m" << std::endl;
            std::cout << "\033[32m" << "============================" << "\033[0m" << std::endl;
        }
    }
    void SetUp() { std::cout << "\033[32m" << "[    CASE  ]" << "\033[0m" << " "; }
    void TearDown() {}
};

#endif