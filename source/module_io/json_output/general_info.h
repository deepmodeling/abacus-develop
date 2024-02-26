#include "module_io/input.h"

//Add json objects to gener_info
namespace Json
{
#ifdef __RAPIDJSON
void gen_general_info(Input *input);
#endif
}