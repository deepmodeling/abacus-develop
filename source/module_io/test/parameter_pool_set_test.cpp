#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "module_io/parameter_pool.h"

/**
 * @brief test func:input_parameters_get() and input_parameters_set() in parameter_pool.h
 */
TEST(InputParameterArgsTest, InputTest)
{

    bool b1 = false;
    ModuleIO::InputParameter inputparaBool;
    inputparaBool.type = BOOL;
    inputparaBool.set((void*)&b1);
    ModuleIO::input_parameters["bool"] = inputparaBool;
    bool b2 = *static_cast<bool*>(ModuleIO::input_parameters["bool"].get());
    std::cout << "bool test:" << b2 << "\n";
    EXPECT_EQ(b2, b1);

    int i1 = 2023;
    ModuleIO::InputParameter inputparaInt;
    inputparaInt.type = INT;
    inputparaInt.set((void*)&i1);
    ModuleIO::input_parameters["int"] = inputparaInt;
    int i2 = *static_cast<int*>(ModuleIO::input_parameters["int"].get());
    std::cout << "int test:" << i2 << "\n";
    EXPECT_EQ(i1, i2);

    SimpleString str = "Hello AISI";
    ModuleIO::InputParameter inputpara;
    inputpara.type = STRING;
    inputpara.set((void*)&str);
    ModuleIO::input_parameters["string"] = inputpara;
    // ModuleIO::input_parameters.insert(std::pair<std::string, ModuleIO::InputParameter>("nupdown", inputpara));
    std::string s = static_cast<SimpleString*>(ModuleIO::input_parameters["string"].get())->c_str();
    std::cout << "string test:" << s << "\n";
    EXPECT_EQ(str.c_str(), s);

    double d1 = 3.1415926;
    ModuleIO::InputParameter inputpara2;
    inputpara2.type = DOUBLE;
    inputpara2.set((void*)&d1);
    ModuleIO::input_parameters["double"] = inputpara2;
    double d2 = *static_cast<double*>(ModuleIO::input_parameters["double"].get());
    std::cout << "double test:" << d2 << "\n";
    EXPECT_EQ(d1, d2);

    SimpleVector<int> vecI;
    for (int i = 0; i < 5; i++)
    {
        vecI.push_back(i);
    }
    ModuleIO::InputParameter inputpara3;
    inputpara3.type = VECTOR_I;
    inputpara3.set((void*)&vecI);
    ModuleIO::input_parameters["vectorINT"] = inputpara3;
    SimpleVector<int> vecI2 = *static_cast<SimpleVector<int>*>(ModuleIO::input_parameters["vectorINT"].get());
    std::cout << "Vector_int test:";
    for (int i = 0; i < vecI2.size(); i++)
    {
        EXPECT_EQ(vecI[i], vecI2[i]);
        std::cout << vecI2[i];
    }
    std::cout << "\n";

    SimpleVector<double> vecD;
    for (int i = 0; i < 5; i++)
    {
        double tmp = i + 0.01;
        vecD.push_back(tmp);
    }
    ModuleIO::InputParameter inputpara4;
    inputpara4.type = VECTOR_D;
    inputpara4.set((void*)&vecD);
    ModuleIO::input_parameters["vectorDOUBLE"] = inputpara4;
    SimpleVector<double> vecD2 = *static_cast<SimpleVector<double>*>(ModuleIO::input_parameters["vectorDOUBLE"].get());
    std::cout << "Vector_double test:";
    for (int i = 0; i < vecD2.size(); i++)
    {
        double value1=vecD[i],value2=vecD2[i];
        EXPECT_EQ(value1, value2);
        std::cout << vecD2[i] << ",";
    }
}