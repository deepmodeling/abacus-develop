#ifndef XCTEST_H
#define XCTEST_H
#include "gtest/gtest.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
class XCTest : public testing::Test
{
  public:
    XCTest ()
    {
        PARAM.input.basis_type = "";
        PARAM.input.cal_force = false;
        PARAM.input.cal_stress = false;
    }
};
#endif