//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_CONVERT_H
#define INPUT_CONVERT_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

using namespace std;

namespace Input_Conv
{
void Convert(void);

// Function `parse_expression` is used to parse input parameters as expressions into vectors
// fn (string): expressions such as "3*1 0 2*0.5 3*0"
// arr (vector): stores parsing results, for example, "3*1 0 2*0.5 1*1.5" can be parsed as [1, 1, 1, 0, 0.5, 0.5, 1.5]
template <typename T> void parse_expression(const std::string &fn, std::vector<T> &arr);
} // namespace Input_Conv

#endif // Input_Convert
