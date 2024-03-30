#ifndef GENERAL_INFO_H
#define GENERAL_INFO_H
#include "module_io/input.h"

/**
* @brief In this part of the code to complete the general_info part of the json tree.
*/
namespace Json
{
#ifdef __RAPIDJSON
void gen_general_info(Input *input);
#endif
}
#endif