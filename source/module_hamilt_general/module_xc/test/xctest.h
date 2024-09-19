#ifndef XCTEST_H
#define XCTEST_H
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
class XCTest: public testing::Test
{
    public:
        XCTest()
        {
            PARAM.input.basis_type = "";
            PARAM.input.cal_force = 0;
            PARAM.input.cal_stress = 0;
        }
};
#endif